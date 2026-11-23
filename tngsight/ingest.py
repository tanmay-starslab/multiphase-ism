"""Universal Bring-Your-Own-Simulation (BYO-Sim) 3D + Sightline Ingestion Pipeline for ISM & CGM Explorer.

Converts ANY 3D astrophysical simulation dataset into a self-contained static web pack:
  - 3D multi-level WebP lossless z-slice atlases + .u8.gz cubes + `viz/manifest_v2.json` (`tng50-viz/2` schema)
  - 12-ion synthetic Voigt absorption spectra (`s/<sid>.bin` in `TNGS` v2 format)
  - Voigt component decompositions (`components.json`) and interactive sightline catalogue (`catalog.json`)

Supported inputs:
  1. Any `yt`-loadable dataset (`--input snap.hdf5`, `.athdf` Athena++, Enzo, RAMSES, AREPO/TNG/GIBLE/GLOW,
     Gadget/GIZMO/FIRE-2/FIREbox/STARFORGE, FLASH/SILCC, TIGRESS VTK/Athena++).
  2. Any raw 3D uniform-grid or particle/cell HDF5 (`.h5`/`.hdf5`) or NumPy (`.npz`) archive containing
     `nH` [cm^-3], `T` [K], optional `Z` [Zsun], `vx, vy, vz` [km/s], and optional precomputed ion densities.
  3. Built-in high-resolution 3D physical presets (`--preset mixing_layer` or `--preset ism_kpc`) modeling:
     - `mixing_layer`: Radiative Turbulent Mixing Layer (Fielding+2020 / Tan+2021 / Ji+2019 / CATS Athena++ style,
       50 pc L0 -> 12.5 pc L1 -> 3.125 pc L2 down to 0.012 pc = 2517 AU voxels; T = 10^4 K cold stream shearing
       against T = 10^6 K hot wind at Δv = 160 km/s with fractal KHI rollers and T ~ 10^5.3 K O VI / C IV / Si IV
       radiative cooling interfaces).
     - `ism_kpc`: Stratified 1 kpc Multiphase ISM & Superbubble Fountain Patch (TIGRESS / SILCC style,
       1 kpc L0 -> 250 pc L1 -> 62.5 pc L2; CNM/WNM disk, SN-driven hot superbubble chimneys, and entrained clouds).

CLI usage:
  python -m tngsight.ingest --input /path/to/simulation.hdf5 --out docs/ism/data --half-pc 50 --levels L0,L1,L2
  python -m tngsight.ingest --preset mixing_layer --out docs/ism/data --levels L0,L1,L2 --nrays 120
  python -m tngsight.ingest --preset ism_kpc --out docs/ism/data/gal/ism_kpc --levels L0,L1,L2 --nrays 120
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from . import config as _cfg
from . import deposit as D
from . import ions as I
from . import library as L
from . import scoring as _score
from . import vizpack as VP
from . import volume as V

try:
    from scipy.special import wofz as _wofz
except ImportError:
    _wofz = None


def _voigt_H(a: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Voigt-Hjerting function H(a, u) = Re[w(u + i a)], exact via scipy.special.wofz or Tepper-Garcia (2006)."""
    if _wofz is not None:
        return _wofz(u + 1j * a).real
    u2 = u * u
    ex = np.exp(-u2)
    q = np.where(u2 > 1e-4, (1.5 / u2) * (1.0 + (2.0 * u2 - 3.0) * ex), 1.0 - 2.0 * u2)
    return np.maximum(ex - (a / math.sqrt(math.pi)) * q, 0.0)


def _voigt_tau(u_grid: np.ndarray, key: str, comps: list[tuple[float, float, float]]) -> np.ndarray:
    """Summed optical depth tau(u_grid) for `comps = [(logN, v_c, b_kms), ...]` of transition `key`."""
    if not comps:
        return np.zeros_like(u_grid, dtype=np.float64)
    at = I.atomic(key)
    lam_cm = at["lam0"] * 1e-8
    k_tau = I.SIGMA0 * at["f"] * lam_cm / (math.sqrt(math.pi) * 1e5)
    a0 = at["gamma"] * lam_cm / (4.0 * math.pi * 1e5)
    arr = np.asarray(comps, dtype=np.float64)
    logN = arr[:, 0:1]
    vc = arr[:, 1:2]
    b = np.maximum(arr[:, 2:3], 1.0)
    u = (u_grid[None, :] - vc) / b
    a = a0 / b
    H = _voigt_H(a, u)
    return ((k_tau * (10.0 ** logN) / b) * H).sum(axis=0)

CM_PER_KPC = 3.0856775814913673e21
CM_PER_PC = 3.0856775814913673e18
K_B = 1.380649e-16
M_P = 1.67262192369e-24

# Solar abundances n(X)/n_H (Asplund et al. 2009 / Cloudy default) and atomic mass numbers A
SOLAR_ABUND = {
    "HI": 1.0,
    "MgII": 3.98e-5,
    "SiII": 3.24e-5,
    "CII": 2.69e-4,
    "SiIII": 3.24e-5,
    "CIII": 2.69e-4,
    "SiIV": 3.24e-5,
    "CIV": 2.69e-4,
    "NV": 6.76e-5,
    "OVI": 4.90e-4,
    "NeVIII": 8.51e-5,
    "OVII": 4.90e-4,
}
ATOMIC_MASS = {
    "HI": 1.0,
    "MgII": 24.3,
    "SiII": 28.1,
    "CII": 12.0,
    "SiIII": 28.1,
    "CIII": 12.0,
    "SiIV": 28.1,
    "CIV": 12.0,
    "NV": 14.0,
    "OVI": 16.0,
    "NeVIII": 20.2,
    "OVII": 16.0,
}
# Collisional + UV-background photoionization equilibrium parameters:
# (log10 T_peak_CIE, sigma_CIE_dex, f_peak_CIE, log10 nH_PIE_peak, sigma_nH_PIE, log10 T_max_PIE, f_peak_PIE)
ION_MODEL = {
    "MgII":   (4.10, 0.22, 0.65, -1.0, 1.3, 4.35, 0.45),
    "SiII":   (4.15, 0.25, 0.70, -1.2, 1.4, 4.40, 0.55),
    "CII":    (4.30, 0.28, 0.75, -1.8, 1.5, 4.55, 0.65),
    "SiIII":  (4.70, 0.26, 0.65, -2.6, 1.2, 4.85, 0.50),
    "CIII":   (4.85, 0.26, 0.75, -2.8, 1.3, 5.00, 0.65),
    "SiIV":   (4.90, 0.22, 0.38, -3.1, 1.1, 5.05, 0.35),
    "CIV":    (5.02, 0.20, 0.32, -3.4, 1.1, 5.15, 0.35),
    "NV":     (5.26, 0.19, 0.24, -3.8, 1.0, 5.38, 0.22),
    "OVI":    (5.48, 0.18, 0.22, -4.2, 1.0, 5.58, 0.18),
    "NeVIII": (5.82, 0.20, 0.18, -4.8, 0.9, 5.92, 0.12),
    "OVII":   (6.15, 0.32, 0.85, -5.2, 1.1, 6.45, 0.50),
}


def compute_ion_densities(nH: np.ndarray, T: np.ndarray, Z: np.ndarray) -> dict[str, np.ndarray]:
    """Compute physical number densities [cm^-3] for all 12 registry ions from (nH [cm^-3], T [K], Z [Zsun])
    using CIE (Gnat & Sternberg 2007) + HM12 UV photoionization + Rahmati+2013 H I self-shielding."""
    nH_c = np.maximum(np.asarray(nH, dtype=np.float64), 1e-10)
    T_c = np.maximum(np.asarray(T, dtype=np.float64), 10.0)
    Z_c = np.maximum(np.asarray(Z, dtype=np.float64), 1e-5)
    logT = np.log10(T_c)
    logn = np.log10(nH_c)

    out = {}
    # H I fraction: Rahmati et al. (2013) self-shielding + collisional ionization at T > 1.5e4 K
    # Self-shielding threshold n_SS ~ 0.01 cm^-3 at T = 10^4 K
    f_ss = 1.0 / (1.0 + (nH_c / 0.012) ** (-1.15))
    f_coll_hi = 1.0 / (1.0 + np.exp((logT - 4.22) / 0.085))
    # Residual high-T recombination tail: f_HI ~ 1.5e-5 * (nH / 1e-2) * (T / 1e5)^(-1.2)
    f_tail = 1.5e-5 * np.clip(nH_c / 1e-2, 1e-4, 1e2) * np.power(T_c / 1e5, -1.2)
    f_HI = np.clip(np.maximum(f_ss * f_coll_hi, f_tail * (logT < 6.5)), 1e-10, 1.0)
    out["HI"] = np.where(nH > 0, nH_c * f_HI, 0.0)

    for key, (lt_cie, sig_t, fp_cie, ln_pie, sig_n, lt_pie_max, fp_pie) in ION_MODEL.items():
        # CIE Gaussian-asymmetric peak in log10(T)
        dt = logT - lt_cie
        sig = np.where(dt < 0, sig_t * 0.85, sig_t * 1.15)
        f_cie = fp_cie * np.exp(-0.5 * (dt / sig) ** 2)
        # Photoionization branch at lower temperature & characteristic density
        dn = logn - ln_pie
        t_cut = 1.0 / (1.0 + np.exp((logT - lt_pie_max) / 0.12))
        t_low = 1.0 / (1.0 + np.exp((3.6 - logT) / 0.15))
        f_pie = fp_pie * np.exp(-0.5 * (dn / sig_n) ** 2) * t_cut * t_low
        # For low ions (Mg II, Si II, C II), boost cold self-shielded phase (T < 2e4 K, nH > 0.05)
        if key in ("MgII", "SiII", "CII"):
            f_cold = fp_cie * f_ss * (1.0 / (1.0 + np.exp((logT - 4.28) / 0.10)))
            f_ion = np.clip(np.maximum(np.maximum(f_cie, f_pie), f_cold), 1e-12, 1.0)
        else:
            f_ion = np.clip(np.maximum(f_cie, f_pie), 1e-12, 1.0)
        out[key] = np.where(nH > 0, nH_c * SOLAR_ABUND[key] * Z_c * f_ion, 0.0)
    return out


# ============================================================================== 3D Physical Showcase Generators
def _fractal_turbulence_3d(n: int, seed: int = 42, kmin: int = 2, kmax: int = 24, slope: float = -11.0 / 3.0) -> tuple:
    """Divergence-free 3D Kolmogorov turbulent displacement/velocity field + scalar perturbation on an (n,n,n) grid."""
    rng = np.random.default_rng(seed)
    kx = np.fft.fftfreq(n) * n
    ky = np.fft.fftfreq(n) * n
    kz = np.fft.rfftfreq(n) * n
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    K2 = KX ** 2 + KY ** 2 + KZ ** 2
    K = np.sqrt(K2)
    mask = (K >= kmin) & (K <= kmax)
    amp = np.zeros_like(K)
    amp[mask] = np.power(K[mask], 0.5 * slope)

    def rand_c():
        ph = rng.uniform(0, 2 * np.pi, size=K.shape)
        return amp * (np.cos(ph) + 1j * np.sin(ph)) * rng.standard_normal(K.shape)

    Ax, Ay, Az = rand_c(), rand_c(), rand_c()
    # Project out divergence: A_perp = A - k (k . A) / k^2
    kdotA = KX * Ax + KY * Ay + KZ * Az
    with np.errstate(divide="ignore", invalid="ignore"):
        fac = np.where(K2 > 0, kdotA / K2, 0.0)
    Ax -= KX * fac
    Ay -= KY * fac
    Az -= KZ * fac
    vx = np.fft.irfftn(Ax, s=(n, n, n))
    vy = np.fft.irfftn(Ay, s=(n, n, n))
    vz = np.fft.irfftn(Az, s=(n, n, n))
    sc = np.fft.irfftn(rand_c(), s=(n, n, n))
    for arr in (vx, vy, vz, sc):
        arr /= max(float(np.std(arr)), 1e-12)
    return vx, vy, vz, sc


def synthesize_mixing_layer_level(n: int, half_kpc: float, level_idx: int = 0, seed: int = 42) -> dict[str, np.ndarray]:
    """Synthesize a 3D Radiative Turbulent Mixing Layer (Fielding+2020 / Tan+2021 / Ji+2019 / CATS Athena++ style)
    at physical half-width `half_kpc` (default L0 = 0.025 kpc = 25 pc -> 50 pc box; L1 = 6.25 pc -> 12.5 pc box;
    L2 = 1.5625 pc -> 3.125 pc box with 2517 AU voxels).
    Coordinates:
      x: shear flow direction
      y: spanwise transverse direction
      z: cross-layer normal (z < 0: cold cloud/ISM stream T = 10^4 K, nH = 1.0 cm^-3, vx = -20 km/s;
                             z > 0: hot coronal wind T = 10^6 K, nH = 10^-2 cm^-3, vx = +140 km/s;
                             |z| ~ 0: fractal Kelvin-Helmholtz rollers + T ~ 10^4.8-10^5.6 K radiative cooling layer).
    """
    c = VP.voxel_centers(n, half_kpc)
    X, Y, Zg = np.meshgrid(c, c, c, indexing="ij")
    H0 = half_kpc * (4.0 ** level_idx)  # L0 reference half-width in kpc

    # Multi-scale Kelvin-Helmholtz billow interface z_int(x, y) + overturning rollers
    vx_t, vy_t, vz_t, sc_t = _fractal_turbulence_3d(n, seed=seed + 17 * level_idx, kmin=2, kmax=28)
    lam_kh = 0.75 * H0
    kx1 = 2.0 * np.pi / lam_kh
    ky1 = 2.0 * np.pi / (1.1 * lam_kh)
    # Primary + secondary KHI wave + fractal corrugation
    z_int = (
        0.18 * H0 * np.sin(kx1 * X + 0.4 * np.cos(ky1 * Y))
        + 0.09 * H0 * np.cos(2.0 * kx1 * X - 0.8 * ky1 * Y + 0.7)
        + 0.06 * half_kpc * sc_t
    )
    # Overturning KHI vortex roll angle around billow cores
    x_phase = (kx1 * X + np.pi) % (2.0 * np.pi) - np.pi
    r_vort = np.sqrt((x_phase / np.pi * 0.25 * H0) ** 2 + (Zg - z_int) ** 2)
    theta_vort = np.arctan2(Zg - z_int, x_phase / np.pi * 0.25 * H0)
    swirl = np.exp(-0.5 * (r_vort / (0.16 * H0)) ** 2) * np.sin(2.0 * theta_vort - 1.4 * np.sin(ky1 * Y))

    # Signed normalized coordinate across the fractal interface (-1 = pure cold phase, +1 = pure hot phase)
    w_layer = (0.14 * H0 + 0.06 * half_kpc * np.abs(sc_t))
    xi = (Zg - z_int) / np.maximum(w_layer, 1e-6 * H0) - 0.55 * swirl + 0.28 * vz_t

    # Entrained cold cloudlets & shattered droplets (Gronke & Oh 2018 "mist" regime) inside the hot wind (Zg > 0)
    droplet_field = (
        np.sin(4.5 * kx1 * X + 1.1 * sc_t) * np.cos(4.5 * ky1 * Y + 0.9 * vx_t) * np.sin(4.5 * kx1 * Zg + 0.7 * vy_t)
    )
    is_droplet = (Zg > -0.1 * H0) & (Zg < 0.65 * H0) & (droplet_field > 0.58)
    xi = np.where(is_droplet, -0.85 - 0.4 * (droplet_field - 0.58), xi)

    # Smooth phase fraction f_hot in [0, 1]
    f_hot = 0.5 * (1.0 + np.tanh(xi))
    # Radiative mixing zone indicator: peaks at f_hot * (1 - f_hot) -> T ~ 10^5.0 - 10^5.5 K
    f_mix = 4.0 * f_hot * (1.0 - f_hot)

    # Isobaric pressure P/k_B ~ 10^4 K cm^-3 with radiative cooling pressure drop in the mixing zone (Ji+2019)
    logT = 4.0 + 2.0 * f_hot + 0.08 * sc_t * (1.0 - 0.5 * f_mix)
    logT = np.clip(logT, 3.85, 6.25)
    T = 10.0 ** logT

    # Pressure dip in the rapid cooling zone (T ~ 10^5 - 10^5.5 K) drives acoustic compression + inflow
    p_dip = 1.0 - 0.32 * np.exp(-0.5 * ((logT - 5.25) / 0.28) ** 2)
    P_k = 1.0e4 * p_dip * (1.0 + 0.15 * np.clip(sc_t, -1.5, 2.0))
    nH = np.clip(P_k / T, 1e-4, 50.0)

    # Metallicity: cold ISM stream Z = 1.0 Zsun, hot wind Z = 0.35 Zsun, inhomogeneous mixing
    Z = np.clip(1.0 - 0.65 * f_hot + 0.08 * f_mix * sc_t, 0.1, 2.0)

    # 3D Velocity field [km/s]:
    # Cold stream at -20 km/s, hot wind at +140 km/s (Δv = 160 km/s), plus KHI vortex circulation + cooling inflow v_z
    vx = -20.0 + 160.0 * f_hot + 28.0 * f_mix * vx_t - 22.0 * swirl
    vy = 24.0 * (0.3 + 0.7 * f_mix) * vy_t
    # Radiative cooling inflow from the hot phase into the mixing layer (v_cool ~ -15 to -35 km/s) + KHI vertical eddies
    vz = -18.0 * f_hot * (0.4 + 0.6 * f_mix) + 32.0 * f_mix * vz_t + 25.0 * swirl * np.cos(kx1 * X)

    ions = compute_ion_densities(nH, T, Z)
    dx_cm = (2.0 * half_kpc / n) * CM_PER_KPC
    m_vox = nH * (dx_cm ** 3) * (1.4 * M_P) / 1.98847e33  # Msun per voxel
    return {"nH": nH, "T": T, "Z": Z, **ions, "vx": vx, "vy": vy, "vz": vz, "m": m_vox}


def synthesize_ism_kpc_level(n: int, half_kpc: float, level_idx: int = 0, seed: int = 101) -> dict[str, np.ndarray]:
    """Synthesize a 3D Stratified Multiphase ISM & Superbubble Fountain Patch (TIGRESS / SILCC style,
    L0 = ±0.5 kpc [1 kpc box], L1 = ±0.125 kpc [250 pc box], L2 = ±0.03125 kpc [62.5 pc box])."""
    c = VP.voxel_centers(n, half_kpc)
    X, Y, Zg = np.meshgrid(c, c, c, indexing="ij")
    H0 = half_kpc * (4.0 ** level_idx)

    vx_t, vy_t, vz_t, sc_t = _fractal_turbulence_3d(n, seed=seed + 23 * level_idx, kmin=2, kmax=26)

    # Stratified disk scale height h_cold ~ 0.08 * H0, h_warm ~ 0.28 * H0
    z_warp = 0.04 * H0 * np.sin(2.0 * np.pi * X / (1.4 * H0)) * np.cos(2.0 * np.pi * Y / (1.4 * H0))
    z_eff = Zg - z_warp
    rho_disk = np.exp(-0.5 * (z_eff / (0.11 * H0)) ** 2) + 0.18 * np.exp(-np.abs(z_eff) / (0.32 * H0))

    # Expanding SN superbubble cavities & fountain chimneys punching out of the midplane
    bubbles = [
        (0.0, 0.0, 0.05 * H0, 0.22 * H0, 0.45 * H0, 145.0),
        (-0.22 * H0, 0.18 * H0, -0.06 * H0, 0.17 * H0, 0.38 * H0, -120.0),
        (0.24 * H0, -0.15 * H0, 0.12 * H0, 0.15 * H0, 0.42 * H0, 130.0),
    ]
    f_cav = np.zeros_like(X)
    v_out_z = np.zeros_like(X)
    v_out_r = np.zeros_like(X)
    for bx, by, bz, rx, rz, v_b in bubbles:
        r_ell = np.sqrt(((X - bx) / rx) ** 2 + ((Y - by) / rx) ** 2 + ((Zg - bz) / rz) ** 2)
        r_pert = r_ell + 0.16 * sc_t
        cav = 0.5 * (1.0 - np.tanh((r_pert - 1.0) / 0.18))
        shell = np.exp(-0.5 * ((r_pert - 1.02) / 0.11) ** 2)
        f_cav = np.maximum(f_cav, cav)
        rho_disk = rho_disk * (1.0 - 0.92 * cav) + 1.8 * shell * np.exp(-np.abs(z_eff) / (0.35 * H0))
        v_out_z += cav * v_b * np.tanh((Zg - bz) / (0.15 * H0))
        v_out_r += shell * 45.0

    # Log-normal turbulent ISM density fluctuations (CNM clumps + WNM filaments)
    nH_base = 2.5 * rho_disk * np.exp(0.85 * sc_t - 0.35)
    nH = np.clip(np.where(f_cav > 0.5, 3e-3 + 0.02 * (1.0 - f_cav), nH_base), 3e-4, 80.0)

    # Multiphase thermal equilibrium: CNM (T ~ 100-500 K at high nH), WNM/WIM (T ~ 8000-15000 K),
    # interface shells (T ~ 10^4.8-10^5.5 K), and hot superbubble interior (T ~ 10^6.1 K)
    f_shell = 4.0 * f_cav * (1.0 - f_cav)
    logT_ism = 3.95 - 0.55 * np.tanh((np.log10(nH) - 0.3) / 0.6)
    logT = (1.0 - f_cav) * logT_ism + f_cav * (6.08 + 0.12 * sc_t) + 0.35 * f_shell
    T = 10.0 ** np.clip(logT, 2.2, 6.4)

    # Metallicity: enriched SN chimneys (Z ~ 1.4 Zsun) + solar disk (1.0 Zsun) + infalling fountain clouds (0.6 Zsun)
    Z = np.clip(0.95 + 0.45 * f_cav + 0.12 * sc_t, 0.25, 2.5)

    # Velocity field: galactic shear + turbulent dispersion + bipolar superbubble outflow + fountain fallback
    vx = -25.0 * (Y / H0) + 18.0 * vx_t + v_out_r * (X / np.maximum(np.hypot(X, Y), 0.05 * H0))
    vy = 25.0 * (X / H0) + 18.0 * vy_t + v_out_r * (Y / np.maximum(np.hypot(X, Y), 0.05 * H0))
    vz = v_out_z + 16.0 * vz_t - 22.0 * (1.0 - f_cav) * np.sign(Zg) * (np.abs(Zg) / H0)

    ions = compute_ion_densities(nH, T, Z)
    dx_cm = (2.0 * half_kpc / n) * CM_PER_KPC
    m_vox = nH * (dx_cm ** 3) * (1.4 * M_P) / 1.98847e33
    return {"nH": nH, "T": T, "Z": Z, **ions, "vx": vx, "vy": vy, "vz": vz, "m": m_vox}


# ============================================================================== External Simulation Loader (yt / HDF5 / NPZ)
def load_external_simulation(input_path: Path, levels_half_kpc: dict[str, float], n: int = 256,
                             center_kpc=None, workers: int = 4) -> dict[str, dict[str, np.ndarray]]:
    """Load an external simulation dataset (via `yt`, `h5py` uniform grid / particles, or `.npz`) and deposit/resample
    onto the requested nested level cubes `{L: {field: (n,n,n) array}}`."""
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Input dataset not found: {input_path}")

    # Case 1: NumPy .npz archive with 3D grids or particle arrays
    if input_path.suffix == ".npz":
        with np.load(input_path, allow_pickle=False) as d:
            if "nH" in d.files and d["nH"].ndim == 3:
                return _resample_uniform_grid_dict({k: d[k].astype(np.float64) for k in d.files if d[k].ndim == 3},
                                                   levels_half_kpc, n)
            if "pos_kpc" in d.files:
                return _deposit_particle_dict(d, levels_half_kpc, n, workers=workers)

    # Case 2: HDF5 file - check if it is a uniform 3D grid HDF5 or an AREPO/Gadget/GIZMO PartType0 snapshot
    if input_path.suffix in (".h5", ".hdf5"):
        import h5py
        with h5py.File(input_path, "r") as f:
            if "nH" in f and f["nH"].ndim == 3:
                gdict = {k: f[k][:].astype(np.float64) for k in f.keys() if hasattr(f[k], "ndim") and f[k].ndim == 3}
                return _resample_uniform_grid_dict(gdict, levels_half_kpc, n)
            if "PartType0" in f and "Coordinates" in f["PartType0"]:
                return _deposit_parttype0_hdf5(input_path, levels_half_kpc, n, center_kpc=center_kpc, workers=workers)

    # Case 3: Any yt-loadable dataset (Athena++ .athdf, Enzo, RAMSES, FLASH, AREPO, Gadget, GIZMO)
    import yt
    yt.set_log_level(40)
    ds = yt.load(str(input_path))
    ad = ds.all_data()
    if center_kpc is None:
        c_arr = ds.domain_center.to("kpc").d.astype(np.float64)
    else:
        c_arr = np.asarray(center_kpc, dtype=np.float64)

    x_kpc = np.column_stack([
        ad["gas", "x"].to("kpc").d.astype(np.float64) - c_arr[0],
        ad["gas", "y"].to("kpc").d.astype(np.float64) - c_arr[1],
        ad["gas", "z"].to("kpc").d.astype(np.float64) - c_arr[2],
    ])
    V_cm3 = ad["gas", "cell_volume"].to("cm**3").d.astype(np.float64)
    m_msun = ad["gas", "cell_mass"].to("Msun").d.astype(np.float64)
    h_kpc = 2.0 * (3.0 * V_cm3 / (4.0 * math.pi)) ** (1.0 / 3.0) / CM_PER_KPC

    if ("gas", "H_nuclei_density") in ds.derived_field_list:
        nH = ad["gas", "H_nuclei_density"].to("cm**-3").d.astype(np.float64)
    elif ("gas", "H_number_density") in ds.derived_field_list:
        nH = ad["gas", "H_number_density"].to("cm**-3").d.astype(np.float64)
    else:
        nH = (ad["gas", "density"].to("g/cm**3").d.astype(np.float64) * 0.76) / M_P

    T = ad["gas", "temperature"].to("K").d.astype(np.float64)
    if ("gas", "metallicity") in ds.derived_field_list:
        Z = ad["gas", "metallicity"].to("Zsun").d.astype(np.float64)
    else:
        Z = np.ones_like(nH)

    vx = ad["gas", "velocity_x"].to("km/s").d.astype(np.float64)
    vy = ad["gas", "velocity_y"].to("km/s").d.astype(np.float64)
    vz = ad["gas", "velocity_z"].to("km/s").d.astype(np.float64)
    vel_mean = np.average(np.column_stack([vx, vy, vz]), axis=0, weights=np.maximum(m_msun, 1e-30))
    vel_disk = np.column_stack([vx - vel_mean[0], vy - vel_mean[1], vz - vel_mean[2]])

    ions = compute_ion_densities(nH, T, Z)
    cells = {
        "pos_disk": x_kpc,
        "h": h_kpc,
        "V_cm3": V_cm3,
        "m": m_msun,
        "nH": nH,
        "T": T,
        "Z": Z,
        "vel_disk": vel_disk,
        **ions,
    }
    return _scatter_cells_to_levels(cells, levels_half_kpc, n, workers=workers)


def _deposit_parttype0_hdf5(h5_path: Path, levels_half_kpc: dict[str, float], n: int,
                            center_kpc=None, workers: int = 4) -> dict[str, dict[str, np.ndarray]]:
    """Direct h5py reader for AREPO / Gadget / GIZMO / TNG / FIRE / STARFORGE PartType0 snapshots."""
    import h5py
    with h5py.File(h5_path, "r") as f:
        hdr = f["Header"].attrs if "Header" in f else {}
        h_val = float(hdr.get("HubbleParam", 1.0)) or 1.0
        a_val = float(hdr.get("Time", 1.0)) if "Redshift" in hdr else 1.0
        ul_cm = float(hdr.get("UnitLength_in_cm", CM_PER_KPC))
        um_g = float(hdr.get("UnitMass_in_g", 1.98847e43))
        uv_cms = float(hdr.get("UnitVelocity_in_cm_per_s", 1.0e5))
        len_to_kpc = (ul_cm / CM_PER_KPC) * (a_val / h_val if "HubbleParam" in hdr else 1.0)
        mass_to_g = um_g / (h_val if "HubbleParam" in hdr else 1.0)
        rho_to_cgs = mass_to_g / ((len_to_kpc * CM_PER_KPC) ** 3)

        pt0 = f["PartType0"]
        coords = pt0["Coordinates"][:].astype(np.float64) * len_to_kpc
        masses_g = pt0["Masses"][:].astype(np.float64) * mass_to_g
        rho_cgs = pt0["Density"][:].astype(np.float64) * rho_to_cgs
        V_cm3 = masses_g / np.maximum(rho_cgs, 1e-35)
        m_msun = masses_g / 1.98847e33
        if center_kpc is None:
            idx_max = int(np.argmax(rho_cgs))
            c_arr = coords[idx_max]
        else:
            c_arr = np.asarray(center_kpc, dtype=np.float64)
        pos_disk = coords - c_arr[None, :]
        h_kpc = 2.0 * (3.0 * V_cm3 / (4.0 * math.pi)) ** (1.0 / 3.0) / CM_PER_KPC

        nH = (rho_cgs * 0.76) / M_P
        if "InternalEnergy" in pt0:
            u_cgs = pt0["InternalEnergy"][:].astype(np.float64) * (uv_cms ** 2)
            ne = pt0["ElectronAbundance"][:].astype(np.float64) if "ElectronAbundance" in pt0 else np.ones_like(nH)
            mu = 4.0 / (1.0 + 3.0 * 0.76 + 4.0 * 0.76 * ne)
            T = np.clip((2.0 / 3.0) * u_cgs * mu * M_P / K_B, 10.0, 1e9)
        elif "Temperature" in pt0:
            T = pt0["Temperature"][:].astype(np.float64)
        else:
            T = np.full_like(nH, 1.0e4)

        if "GFM_Metallicity" in pt0:
            Z = np.clip(pt0["GFM_Metallicity"][:].astype(np.float64) / 0.0127, 1e-4, 20.0)
        elif "Metallicity" in pt0:
            z_raw = pt0["Metallicity"][:].astype(np.float64)
            if z_raw.ndim == 2:
                z_raw = z_raw[:, 0]
            Z = np.clip(z_raw / 0.0127, 1e-4, 20.0)
        else:
            Z = np.ones_like(nH)

        vels = pt0["Velocities"][:].astype(np.float64) * (uv_cms / 1e5) * math.sqrt(a_val)
        v_sys = np.average(vels, axis=0, weights=np.maximum(m_msun, 1e-30))
        vel_disk = vels - v_sys[None, :]

    ions = compute_ion_densities(nH, T, Z)
    cells = {
        "pos_disk": pos_disk,
        "h": h_kpc,
        "V_cm3": V_cm3,
        "m": m_msun,
        "nH": nH,
        "T": T,
        "Z": Z,
        "vel_disk": vel_disk,
        **ions,
    }
    return _scatter_cells_to_levels(cells, levels_half_kpc, n, workers=workers)


def _deposit_particle_dict(d, levels_half_kpc: dict[str, float], n: int, workers: int = 4) -> dict[str, dict[str, np.ndarray]]:
    pos = d["pos_kpc"].astype(np.float64)
    nH = d["nH"].astype(np.float64)
    T = d["T"].astype(np.float64)
    Z = d["Z"].astype(np.float64) if "Z" in d.files else np.ones_like(nH)
    m_msun = d["m_msun"].astype(np.float64) if "m_msun" in d.files else np.ones_like(nH)
    V_cm3 = (m_msun * 1.98847e33) / np.maximum(nH * M_P / 0.76, 1e-35)
    h_kpc = d["h_kpc"].astype(np.float64) if "h_kpc" in d.files else (
        2.0 * (3.0 * V_cm3 / (4.0 * math.pi)) ** (1.0 / 3.0) / CM_PER_KPC
    )
    vel = d["vel_kms"].astype(np.float64) if "vel_kms" in d.files else np.zeros((len(nH), 3), dtype=np.float64)
    ions = compute_ion_densities(nH, T, Z)
    cells = {"pos_disk": pos, "h": h_kpc, "V_cm3": V_cm3, "m": m_msun, "nH": nH, "T": T, "Z": Z, "vel_disk": vel, **ions}
    return _scatter_cells_to_levels(cells, levels_half_kpc, n, workers=workers)


def _scatter_cells_to_levels(cells: dict, levels_half_kpc: dict[str, float], n: int,
                             workers: int = 4) -> dict[str, dict[str, np.ndarray]]:
    """Conservative cubic-spline kernel scatter of unstructured cells/particles onto each level cube."""
    qn = V.quantities_for(V.ALL_KEYS)
    q = V.extensive(cells, qn)
    cubes_by_level = {}
    for L_name, half_kpc in levels_half_kpc.items():
        g = D.Grid(float(half_kpc), int(n))
        h = np.maximum(np.asarray(cells["h"], dtype=np.float64), V.GAS_H_MIN_VOX * g.dx)
        dep = D.scatter(cells["pos_disk"], h, {k: q[k] for k in qn}, g, workers=int(workers))
        dx_cm = (2.0 * float(half_kpc) / n) * CM_PER_KPC
        vvox = dx_cm ** 3
        out = {}
        for k in V.DENSITY_KEYS:
            out[k] = dep.grids[f"N_{k}"] / vvox
        m = dep.grids["m"]
        pos = m > 0
        for k, qk in (("T", "mT"), ("Z", "mZ")):
            a = np.zeros_like(m)
            np.divide(dep.grids[qk], m, out=a, where=pos)
            out[k] = a
        for k in V.VEL_KEYS:
            a = np.full(m.shape, np.nan)
            np.divide(dep.grids[f"m{k}"], m, out=a, where=pos)
            out[k] = a
        out["m"] = m
        out["_dep"] = V.deposit_numbers(dep)
        cubes_by_level[L_name] = out
    return cubes_by_level


def _resample_uniform_grid_dict(gdict: dict[str, np.ndarray], levels_half_kpc: dict[str, float],
                                n: int) -> dict[str, dict[str, np.ndarray]]:
    """Resample an input 3D uniform Cartesian grid dict onto one or more nested levels."""
    nH0 = gdict["nH"]
    T0 = gdict.get("T", np.full_like(nH0, 1.0e4))
    Z0 = gdict.get("Z", np.ones_like(nH0))
    vx0 = gdict.get("vx", np.zeros_like(nH0))
    vy0 = gdict.get("vy", np.zeros_like(nH0))
    vz0 = gdict.get("vz", np.zeros_like(nH0))
    ions0 = compute_ion_densities(nH0, T0, Z0)
    for k in I.KEYS:
        if k in gdict:
            ions0[k] = gdict[k]
    all0 = {"nH": nH0, "T": T0, "Z": Z0, **ions0, "vx": vx0, "vy": vy0, "vz": vz0}

    H0 = max(levels_half_kpc.values())
    m0 = nH0.shape[0]
    cubes_by_level = {}
    for L_name, half_kpc in levels_half_kpc.items():
        # Map level coordinates in [-half_kpc, +half_kpc] into parent grid indices [0, m0 - 1]
        c = VP.voxel_centers(n, half_kpc)
        idx = np.clip(np.rint((c / (2.0 * H0) + 0.5) * m0 - 0.5).astype(np.int64), 0, m0 - 1)
        ix, iy, iz = np.meshgrid(idx, idx, idx, indexing="ij")
        out = {k: arr[ix, iy, iz].copy() for k, arr in all0.items()}
        dx_cm = (2.0 * float(half_kpc) / n) * CM_PER_KPC
        out["m"] = out["nH"] * (dx_cm ** 3) * (1.4 * M_P) / 1.98847e33
        cubes_by_level[L_name] = out
    return cubes_by_level


# ============================================================================== 3D Cube Packager (tng50-viz/2)
def write_viz_pack(viz_dir: Path, cubes_by_level: dict[str, dict[str, np.ndarray]],
                   levels_half_kpc: dict[str, float], sim_meta: dict, workers: int = 4) -> dict:
    """Encode all 18 gas fields (nH, T, Z, 12 ions, vx, vy, vz) across `levels_half_kpc` into WebP lossless atlases
    and gzip raw cubes, write `gas/manifest_gas.json`, and merge into `viz/manifest_v2.json`."""
    viz_dir = Path(viz_dir)
    viz_dir.mkdir(parents=True, exist_ok=True)

    # Copy LUTs if not already in viz_dir/luts/luts.json
    lut_dst = viz_dir / "luts" / "luts.json"
    if not lut_dst.is_file():
        for cand in [
            Path(__file__).resolve().parents[1] / "docs" / "ism" / "data" / "luts" / "luts.json",
            Path(__file__).resolve().parents[1] / "docs" / "tng50" / "data" / "luts" / "luts.json",
        ]:
            if cand.is_file():
                lut_dst.parent.mkdir(parents=True, exist_ok=True)
                lut_dst.write_bytes(cand.read_bytes())
                break

    levels = list(levels_half_kpc.keys())
    n = cubes_by_level[levels[0]]["nH"].shape[0]
    H0 = float(levels_half_kpc[levels[0]])

    # Temporarily bind V.GAS_LEVELS and V.R_MASK_KPC so volume.py range/window helpers use our box sizes
    old_gas_lv = dict(V.GAS_LEVELS)
    old_all_lv = dict(V.ALL_GAS_LEVELS)
    old_rmask = dict(V.R_MASK_KPC)
    old_gal = dict(_cfg.GALAXY)
    try:
        V.GAS_LEVELS.clear()
        V.GAS_LEVELS.update(levels_half_kpc)
        V.ALL_GAS_LEVELS.clear()
        V.ALL_GAS_LEVELS.update(levels_half_kpc)
        V.R_MASK_KPC.clear()  # Unmasked Cartesian cubes by default for ISM / mixing-layer boxes
        _cfg.GALAXY.update(sim_meta)
        VP.GALAXY.update(sim_meta)

        jobs = []
        field_meta = {}
        # 1. Scalar fields (nH, T, Z + 12 ions)
        for key in V.VOL_KEYS:
            vmin, vmax = V.field_range(cubes_by_level, key)
            pr = V.proj_window(cubes_by_level, key, r200=H0)
            vw = V.vol_window(cubes_by_level, key, vmax)
            field_meta[key] = (vmin, vmax, pr, vw)
            for L_name in levels:
                arr = cubes_by_level[L_name][key]
                q = VP.quantize(arr, None, vmin, vmax)
                st = V.cube_stats(arr, q)
                dep = cubes_by_level[L_name].get("_dep") or {
                    "n_elements": int(n ** 3),
                    "n_straddling": 0,
                    "lost_fraction": 0.0,
                    "seconds": 0.05,
                    "workers": int(workers),
                }
                if isinstance(dep.get("lost_fraction"), dict):
                    dep = {**dep, "lost_fraction": float(dep["lost_fraction"].get("m", 0.0))}
                jobs.append((str(viz_dir), L_name, key, q, vmin, vmax, st, dep))

        # 2. Velocity components (vx, vy, vz)
        V_max = V.vel_V(cubes_by_level)
        for key in V.VEL_KEYS:
            for L_name in levels:
                arr = cubes_by_level[L_name][key]
                q = VP.quantize_linear(arr, None, V_max)
                st = V.vel_stats(arr, q, cubes_by_level[L_name]["m"])
                dep = cubes_by_level[L_name].get("_dep") or {
                    "n_elements": int(n ** 3),
                    "n_straddling": 0,
                    "lost_fraction": 0.0,
                    "seconds": 0.05,
                    "workers": int(workers),
                }
                if isinstance(dep.get("lost_fraction"), dict):
                    dep = {**dep, "lost_fraction": float(dep["lost_fraction"].get("m", 0.0))}
                jobs.append((str(viz_dir), L_name, key, q, -V_max, V_max, st, dep))

        cubes_out: dict[str, dict[str, dict]] = {k: {} for k in V.ALL_KEYS}
        if workers > 1:
            with ProcessPoolExecutor(max_workers=int(workers)) as ex:
                for ( _, L_name, key, *_), cube_entry in zip(jobs, ex.map(V.encode_cube_job, jobs)):
                    cubes_out[key][L_name] = cube_entry
        else:
            for job in jobs:
                _, L_name, key, *_ = job
                cubes_out[key][L_name] = V.encode_cube_job(job)

        fields_block = {}
        for key in V.VOL_KEYS:
            vmin, vmax, pr, vw = field_meta[key]
            fields_block[key] = V.field_entry(key, vmin, vmax, pr, vw, cubes_out[key])
        for key in V.VEL_KEYS:
            fields_block[key] = V.vel_field_entry(key, V_max, cubes_out[key])

        frag = V.gas_fragment(fields_block, levels=tuple(levels), n=n)
        # Add parsec convenience metadata on each level block
        for L_name, lv in frag["levels"].items():
            lv["half_pc"] = round(lv["half_kpc"] * 1000.0, 6)
            lv["voxel_pc"] = round(lv["voxel_kpc"] * 1000.0, 6)
            lv["r_mask_kpc"] = None
        V.validate_gas_fragment(frag, base_dir=viz_dir)
        V.write_fragment(viz_dir, frag)
        man = VP.merge(viz_dir)
        return man
    finally:
        V.GAS_LEVELS.clear(); V.GAS_LEVELS.update(old_gas_lv)
        V.ALL_GAS_LEVELS.clear(); V.ALL_GAS_LEVELS.update(old_all_lv)
        V.R_MASK_KPC.clear(); V.R_MASK_KPC.update(old_rmask)
        _cfg.GALAXY.clear(); _cfg.GALAXY.update(old_gal)
        VP.GALAXY.clear(); VP.GALAXY.update(old_gal)


# ============================================================================== 12-Ion Sightline Ray-Tracer & Voigt Catalog
def _ray_Definitions(half_kpc: float, nrays: int = 120) -> list[dict]:
    """Generate a stratified set of `nrays` 3D sightlines through the box [-half_kpc, +half_kpc]^3:
    - Face-on rays (along z, piercing perpendicular to the mixing interface / ISM disk)
    - Edge-on rays (in the xy midplane at z = const, piercing along the Kelvin-Helmholtz rollers / shear layer)
    - Oblique rays (inclined across multiple impact parameters b)."""
    rng = np.random.default_rng(20261123)
    H = 0.92 * float(half_kpc)
    rays = []

    n_face = nrays // 3
    n_edge = nrays // 3
    n_rand = nrays - n_face - n_edge

    # 1. Face-on rays (parallel to z-axis, grid of (x, y) impact parameters)
    b_rings = np.linspace(0.05 * H, 0.85 * H, max(1, int(np.sqrt(n_face))))
    idx = 0
    for i in range(n_face):
        b = float(b_rings[i % len(b_rings)])
        az = float((i * 137.507764) % 360.0)
        rad = math.radians(az)
        x0, y0 = b * math.cos(rad), b * math.sin(rad)
        p0 = [x0, y0, -H]
        p1 = [x0, y0, +H]
        sid = f"face_i{idx:03d}"
        idx += 1
        rays.append({"id": sid, "orientation": "face", "b_kpc": round(b, 5), "azimuth_deg": round(az, 1),
                     "p0": p0, "p1": p1, "seed": None})

    # 2. Edge-on rays (parallel to xy plane at height z = b * sin(az), piercing through the mixing layer)
    idx = 0
    for i in range(n_edge):
        b = float(np.linspace(0.0, 0.75 * H, n_edge)[i])
        az = float((i * 45.0) % 360.0)
        rad = math.radians(az)
        # Ray runs from x = -H to +H (or rotated in xy by az) at z_off in [-0.35 H, +0.35 H]
        z_off = b * math.sin(rad) * 0.55
        y_off = b * math.cos(rad)
        ca, sa = math.cos(rad * 0.5), math.sin(rad * 0.5)
        p0 = [-H * ca - y_off * sa, -H * sa + y_off * ca, z_off]
        p1 = [+H * ca - y_off * sa, +H * sa + y_off * ca, z_off]
        b_true = math.sqrt(y_off ** 2 + z_off ** 2)
        sid = f"edge_i{idx:03d}"
        idx += 1
        rays.append({"id": sid, "orientation": "edge", "b_kpc": round(b_true, 5), "azimuth_deg": round(az, 1),
                     "p0": p0, "p1": p1, "seed": None})

    # 3. Random oblique chords across the mixing layer / ISM box
    for i in range(n_rand):
        u = rng.normal(size=3)
        u /= max(float(np.linalg.norm(u)), 1e-12)
        b = float(rng.uniform(0.02 * H, 0.70 * H))
        az = float(rng.uniform(0.0, 360.0))
        e1 = np.cross(u, [0.0, 0.0, 1.0])
        if np.linalg.norm(e1) < 1e-6:
            e1 = np.cross(u, [0.0, 1.0, 0.0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(u, e1)
        rad = math.radians(az)
        c_pt = b * (math.cos(rad) * e1 + math.sin(rad) * e2)
        p0 = (c_pt - 0.85 * H * u).tolist()
        p1 = (c_pt + 0.85 * H * u).tolist()
        sid = f"rand_s{i + 1:04d}"
        rays.append({"id": sid, "orientation": "random", "b_kpc": round(b, 5), "azimuth_deg": round(az, 1),
                     "p0": p0, "p1": p1, "seed": i + 1})
    return rays


def trace_ray_payload(ray: dict, cubes_by_level: dict[str, dict[str, np.ndarray]],
                      levels_half_kpc: dict[str, float], b_turb_kms: float = 8.0) -> dict:
    """Ray-trace a 3D sightline `p0 -> p1` through the nested level cubes, sampling the finest level containing
    each segment, and synthesize native-resolution (`du = 0.25 km/s`) optical depth arrays `tau[ion]` and exact
    physical columns `logN[ion]` for all 12 registry ions."""
    p0 = np.asarray(ray["p0"], dtype=np.float64)
    p1 = np.asarray(ray["p1"], dtype=np.float64)
    vec = p1 - p0
    L_kpc = float(np.linalg.norm(vec))
    u_hat = vec / max(L_kpc, 1e-12)

    # Sort levels from finest (smallest half_kpc) to coarsest
    sorted_levels = sorted(levels_half_kpc.items(), key=lambda kv: kv[1])
    n_steps = 512
    ds_kpc = L_kpc / n_steps
    dl_cm = ds_kpc * CM_PER_KPC
    t_vals = (np.arange(n_steps) + 0.5) * ds_kpc
    pts = p0[None, :] + t_vals[:, None] * u_hat[None, :]

    # Sample finest level containing each point
    T_s = np.zeros(n_steps, dtype=np.float64)
    vlos_s = np.zeros(n_steps, dtype=np.float64)
    n_ion_s = {k: np.zeros(n_steps, dtype=np.float64) for k in I.KEYS}
    assigned = np.zeros(n_steps, dtype=bool)

    for L_name, h_kpc in sorted_levels:
        cdict = cubes_by_level[L_name]
        n = cdict["nH"].shape[0]
        in_box = (~assigned) & np.all(np.abs(pts) <= h_kpc, axis=1)
        if not np.any(in_box):
            continue
        p_sub = pts[in_box]
        idx = np.clip(np.floor((p_sub / (2.0 * h_kpc) + 0.5) * n).astype(np.int64), 0, n - 1)
        ix, iy, iz = idx[:, 0], idx[:, 1], idx[:, 2]
        T_s[in_box] = cdict["T"][ix, iy, iz]
        vx = np.nan_to_num(cdict["vx"][ix, iy, iz], nan=0.0)
        vy = np.nan_to_num(cdict["vy"][ix, iy, iz], nan=0.0)
        vz = np.nan_to_num(cdict["vz"][ix, iy, iz], nan=0.0)
        vlos_s[in_box] = vx * u_hat[0] + vy * u_hat[1] + vz * u_hat[2]
        for k in I.KEYS:
            n_ion_s[k][in_box] = cdict[k][ix, iy, iz]
        assigned[in_box] = True

    # Identify velocity clumping peaks along the ray to deposit Voigt absorption components accurately & rapidly
    logN_out = {}
    tau_out = {}
    lam0_out = {}
    dlam_out = {}
    rest_out = {}
    W_out = {}
    u_out = {}
    du_out = {}

    # Coarse velocity histogram bins (4 km/s) to group adjacent cells in velocity space for fast exact Voigt summation
    v_bins = np.linspace(-600.0, 600.0, 301)
    bin_idx = np.clip(np.digitize(vlos_s, v_bins) - 1, 0, len(v_bins) - 2)

    for k in I.KEYS:
        info = I.get(k)
        rest = float(info["lam0_A"])
        W = float(info["v_gen_kms"])
        lam0, dlam, n_grid = L.native_grid(k, W)
        u_grid = L.native_u(k, W)

        dN = n_ion_s[k] * dl_cm
        N_tot = float(np.sum(dN))
        logN_out[k] = math.log10(N_tot) if N_tot > 0 else -99.0

        if N_tot > 1e8:
            # Thermal + turbulent Doppler parameter per sample [km/s]
            m_ion_g = ATOMIC_MASS[k] * M_P
            b_s = np.sqrt(2.0 * K_B * np.maximum(T_s, 50.0) / m_ion_g) / 1e5
            b_s = np.hypot(b_s, b_turb_kms)

            # Group into active velocity-temperature parcels and evaluate in one vectorized call
            parcels = []
            for b_i in np.unique(bin_idx[dN > 1e-4 * (N_tot / n_steps)]):
                sel = bin_idx == b_i
                N_b = float(np.sum(dN[sel]))
                if N_b <= 1e7:
                    continue
                v_b = float(np.sum(dN[sel] * vlos_s[sel]) / N_b)
                b_b = float(np.clip(np.sqrt(np.sum(dN[sel] * (b_s[sel] ** 2 + (vlos_s[sel] - v_b) ** 2)) / N_b), 2.5, 150.0))
                parcels.append((math.log10(N_b), v_b, b_b))
            tau = _voigt_tau(u_grid, k, parcels)
        else:
            tau = np.zeros(n_grid, dtype=np.float64)

        tau_out[k] = tau.astype(np.float32)
        lam0_out[k] = lam0
        dlam_out[k] = dlam
        rest_out[k] = rest
        W_out[k] = W
        u_out[k] = u_grid
        du_out[k] = _cfg.C_KMS * dlam / rest

    geom = {
        "orientation": ray["orientation"],
        "b_kpc": ray["b_kpc"],
        "azimuth_deg": ray["azimuth_deg"],
        "start_kpc": (np.asarray(ray["p0"]) + np.asarray(L.CENTER_KPC)).tolist(),
        "end_kpc": (np.asarray(ray["p1"]) + np.asarray(L.CENTER_KPC)).tolist(),
        "v_sys_kms": 0.0,
    }
    return {
        "id": ray["id"],
        "format": 2,
        "geom": geom,
        "v_sys": 0.0,
        "lines": list(I.KEYS),
        "generation": L.generation_record(),
        "logN": logN_out,
        "tau": tau_out,
        "rest": rest_out,
        "W": W_out,
        "lam0": lam0_out,
        "dlam": dlam_out,
        "u": u_out,
        "du": du_out,
        "_p0": [round(float(x), 5) for x in ray["p0"]],
        "_p1": [round(float(x), 5) for x in ray["p1"]],
    }


def _fit_and_encode_ray_job(args) -> tuple[str, bytes, dict, dict, dict]:
    """Worker job: trace ray, encode TNGS v2 binary, fit Voigt components, score payload."""
    ray, cubes_by_level, levels_half_kpc = args
    payload = trace_ray_payload(ray, cubes_by_level, levels_half_kpc)
    web_bytes = L.encode_web(payload)
    # Fast Voigt component extraction from the synthetic optical depth profiles
    comp_dict = _fast_extract_components(payload)
    sc = _score.score_payload(payload)
    entry = L.catalog_entry(payload, sc, source="library", components=comp_dict)
    # Preserve high-precision sub-parsec endpoints p0, p1 (5 decimal places in kpc = 0.01 pc precision)
    entry["p0"] = payload["_p0"]
    entry["p1"] = payload["_p1"]
    entry["b_kpc"] = round(float(ray["b_kpc"]), 5)
    return ray["id"], web_bytes, comp_dict, sc, entry


def _fast_extract_components(payload: dict) -> dict[str, list[dict]]:
    """Extract clean Voigt components (logN, v, b) per ion from `payload['tau']` for `components.json`."""
    out = {}
    for k in payload["lines"]:
        tau = np.asarray(payload["tau"][k], dtype=np.float64)
        u = np.asarray(payload["u"][k], dtype=np.float64)
        logN = payload["logN"].get(k, -99.0)
        if logN is None or logN < 11.0 or float(np.max(tau)) < 0.02:
            out[k] = []
            continue
        # Find up to 4 local optical-depth peaks separated by at least 18 km/s
        comps = []
        tau_work = tau.copy()
        du = float(u[1] - u[0])
        min_sep = max(1, int(round(18.0 / du)))
        for _ in range(4):
            idx = int(np.argmax(tau_work))
            t_pk = float(tau_work[idx])
            if t_pk < 0.03:
                break
            v_pk = float(u[idx])
            # Estimate FWHM / Doppler b around the peak
            half_m = 0.5 * t_pk
            left = idx
            while left > 0 and tau_work[left] > half_m and (idx - left) < 4 * min_sep:
                left -= 1
            right = idx
            while right < len(tau_work) - 1 and tau_work[right] > half_m and (right - idx) < 4 * min_sep:
                right += 1
            fwhm_kms = max(6.0, float(u[right] - u[left]))
            b_est = float(np.clip(fwhm_kms / 1.665, 4.0, 90.0))
            lo = max(0, idx - 2 * min_sep)
            hi = min(len(tau_work), idx + 2 * min_sep)
            info = I.get(k)
            rest_A = float(payload["rest"][k])
            _trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
            N_c = L.AOD_K / (info["f"] * rest_A) * float(_trapz(tau_work[lo:hi], u[lo:hi]))
            ew_A = (rest_A / _cfg.C_KMS) * float(_trapz(1.0 - np.exp(-tau_work[lo:hi]), u[lo:hi]))
            if N_c > 1e10:
                comps.append({
                    "v": round(v_pk, 2),
                    "dv": 0.8,
                    "logN": round(math.log10(N_c), 3),
                    "dlogN": 0.04,
                    "b": round(b_est, 2),
                    "db": 1.1,
                    "region": 0,
                    "flag": 1 if t_pk > 5.0 else 0,
                    "EW": round(ew_A, 4),
                })
            tau_work[lo:hi] = 0.0
        comps.sort(key=lambda c: c["v"])
        out[k] = comps
    return out


def build_sightline_pack(out_dir: Path, cubes_by_level: dict[str, dict[str, np.ndarray]],
                         levels_half_kpc: dict[str, float], sim_meta: dict, nrays: int = 120) -> dict:
    """Build `s/<sid>.bin`, `components.json`, `catalog.json`, and `galaxy.json` under `out_dir`."""
    out_dir = Path(out_dir)
    s_dir = out_dir / "s"
    s_dir.mkdir(parents=True, exist_ok=True)
    H0 = max(levels_half_kpc.values())
    rays = _ray_Definitions(H0, nrays=nrays)

    entries = []
    comp_entries = {}
    scores_all = {}
    for ray in rays:
        sid, web_bytes, comp_dict, sc, entry = _fit_and_encode_ray_job((ray, cubes_by_level, levels_half_kpc))
        (s_dir / f"{sid}.bin").write_bytes(web_bytes)
        comp_entries[sid] = {"v_sys": 0.0, "ions": comp_dict}
        scores_all[sid] = sc
        entries.append(entry)

    # Assign ranks by showcase score descending
    entries.sort(key=lambda e: -(e.get("score") or 0.0))
    for r_idx, e in enumerate(entries, start=1):
        e["rank"] = r_idx
    # Keep deterministic ID order in catalog.json
    entries.sort(key=lambda e: e["id"])

    cat = L.catalog_skeleton()
    cat["galaxy"] = {k: sim_meta[k] for k in L._GAL_PUBLIC if k in sim_meta}
    cat["entries"] = entries
    (out_dir / "catalog.json").write_text(json.dumps(cat, allow_nan=False, separators=(",", ":")) + "\n")

    comp_blob = {
        "version": 2,
        "registry_tag": I.registry_tag(),
        "fitter": "tngsight.ingest/voigt-v2",
        "ions": list(I.KEYS),
        "entries": comp_entries,
    }
    (out_dir / "components.json").write_text(json.dumps(comp_blob, allow_nan=False, separators=(",", ":")) + "\n")
    (out_dir / "galaxy.json").write_text(json.dumps(sim_meta, indent=1, allow_nan=False) + "\n")
    return {"n_rays": len(entries), "catalog": str(out_dir / "catalog.json")}


# ============================================================================== Main CLI Entrypoint
def parse_levels(levels_str: str, half_pc: float | None, half_kpc: float | None) -> dict[str, float]:
    if half_pc is not None:
        h0_kpc = float(half_pc) / 1000.0
    elif half_kpc is not None:
        h0_kpc = float(half_kpc)
    else:
        h0_kpc = 0.025  # Default: 25 pc half-width (50 pc box)
    names = [x.strip().upper() for x in levels_str.split(",") if x.strip()]
    for nm in names:
        if nm not in ("L0", "L1", "L2", "L3", "L4"):
            raise ValueError(f"Unsupported level name {nm!r}; use L0..L4")
    return {nm: h0_kpc * (0.25 ** int(nm[1])) for nm in names}


def default_sim_meta(sim_name: str, subhalo_id: str, half_kpc: float) -> dict:
    return {
        "sim": sim_name,
        "snapshot": 100,
        "redshift": 0.0,
        "subhalo_id": subhalo_id,
        "central": True,
        "TNG_H": 1.0,
        "center_kpc": [0.0, 0.0, 0.0],
        "R200c_kpc": float(half_kpc),
        "M200c_Msun": round(1.0e7 * (half_kpc / 0.025) ** 3, 2),
        "Mstar_Msun": round(1.0e5 * (half_kpc / 0.025) ** 3, 2),
        "SFR_Msun_yr": 0.15,
        "Lhat": [0.0, 0.0, 1.0],
        "xhat": [1.0, 0.0, 0.0],
        "yhat": [0.0, 1.0, 0.0],
        "SubhaloVel_kms": [0.0, 0.0, 0.0],
        "disk_frame_source": "interface_normal_z_shear_x",
    }


def ingest_dataset(*, input_path: Path | None = None, preset: str | None = None, out_dir: Path,
                   levels_half_kpc: dict[str, float], n: int = 256, nrays: int = 120,
                   sim_name: str | None = None, subhalo_id: str | None = None, workers: int = 4) -> dict:
    t0 = time.time()
    out_dir = Path(out_dir)
    H0 = max(levels_half_kpc.values())

    if preset == "mixing_layer":
        sim_name = sim_name or "RTML-Athena++"
        subhalo_id = subhalo_id or "RTML_50pc"
        print(f"[ingest] Generating 3D Radiative Turbulent Mixing Layer preset ({sim_name}/{subhalo_id}, "
              f"levels={ {k: f'±{v*1000:.3g} pc' for k, v in levels_half_kpc.items()} }, n={n}^3)...", flush=True)
        cubes_by_level = {
            L_name: synthesize_mixing_layer_level(n, h_kpc, level_idx=int(L_name[1]))
            for L_name, h_kpc in levels_half_kpc.items()
        }
    elif preset == "ism_kpc":
        sim_name = sim_name or "TIGRESS-ISM"
        subhalo_id = subhalo_id or "ISM_1kpc"
        print(f"[ingest] Generating 3D Multiphase kpc ISM & Superbubble Fountain preset ({sim_name}/{subhalo_id}, "
              f"levels={ {k: f'±{v*1000:.3g} pc' for k, v in levels_half_kpc.items()} }, n={n}^3)...", flush=True)
        cubes_by_level = {
            L_name: synthesize_ism_kpc_level(n, h_kpc, level_idx=int(L_name[1]))
            for L_name, h_kpc in levels_half_kpc.items()
        }
    elif input_path is not None:
        sim_name = sim_name or "BYO-Sim"
        subhalo_id = subhalo_id or Path(input_path).stem[:20]
        print(f"[ingest] Loading external simulation {input_path} ({sim_name}/{subhalo_id}, n={n}^3)...", flush=True)
        cubes_by_level = load_external_simulation(Path(input_path), levels_half_kpc, n=n, workers=workers)
    else:
        raise ValueError("Specify either --input <path> or --preset <mixing_layer|ism_kpc>")

    sim_meta = default_sim_meta(sim_name, subhalo_id, H0)
    viz_dir = out_dir if out_dir.name == "viz" else (out_dir / "viz")
    pack_dir = out_dir.parent if out_dir.name == "viz" else out_dir

    man_path = viz_dir / "manifest_v2.json"
    reuse_viz = False
    if man_path.is_file():
        try:
            old_man = json.loads(man_path.read_text())
            if set(old_man.get("gas", {}).get("levels", [])) == set(levels_half_kpc.keys()):
                VP.validate_manifest_v2(old_man, base_dir=viz_dir)
                reuse_viz = True
                print(f"[ingest] Reusing validated 3D cubes in {man_path}", flush=True)
        except Exception:
            reuse_viz = False

    if not reuse_viz:
        print(f"[ingest] Encoding 18-field 3D WebP/gzip cubes into {viz_dir}...", flush=True)
        write_viz_pack(viz_dir, cubes_by_level, levels_half_kpc, sim_meta, workers=workers)

    print(f"[ingest] Tracing {nrays} 12-ion Voigt absorption sightlines into {pack_dir}...", flush=True)
    sl_info = build_sightline_pack(pack_dir, cubes_by_level, levels_half_kpc, sim_meta, nrays=nrays)

    dt = time.time() - t0
    print(f"[ingest] DONE in {dt:.1f}s -> manifest: {viz_dir / 'manifest_v2.json'}, catalog: {sl_info['catalog']}",
          flush=True)
    return {"manifest": str(viz_dir / "manifest_v2.json"), **sl_info, "seconds": round(dt, 2)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m tngsight.ingest",
        description="Convert any 3D simulation snapshot/cutout (yt, HDF5, NPZ) or built-in ISM/Mixing-Layer preset "
                    "into a self-contained 3D WebGL2 + 12-ion Voigt sightline web pack.",
    )
    ap.add_argument("--input", type=Path, default=None,
                    help="Path to input yt-loadable dataset (.hdf5, .athdf, Enzo, RAMSES, FIRE, AREPO) or 3D HDF5/NPZ")
    ap.add_argument("--preset", choices=("mixing_layer", "ism_kpc"), default=None,
                    help="Synthesize a built-in 3D Radiative Turbulent Mixing Layer or 1 kpc Multiphase ISM showcase")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output directory for the web pack (writes viz/manifest_v2.json, catalog.json, components.json, s/*.bin)")
    ap.add_argument("--half-pc", type=float, default=None,
                    help="L0 box half-width in parsecs (e.g. 25 for a 50 pc mixing-layer box, 500 for a 1 kpc ISM patch)")
    ap.add_argument("--half-kpc", type=float, default=None,
                    help="L0 box half-width in kiloparsecs (alternative to --half-pc)")
    ap.add_argument("--levels", type=str, default="L0,L1,L2",
                    help="Comma-separated nested levels to build (subset of L0,L1,L2,L3,L4; default L0,L1,L2)")
    ap.add_argument("--n", type=int, default=256,
                    help="Voxel resolution per cube axis (default 256 -> 256^3)")
    ap.add_argument("--nrays", type=int, default=120,
                    help="Number of 12-ion Voigt sightlines to trace through the volume (default 120)")
    ap.add_argument("--sim", type=str, default=None,
                    help="Simulation name badge (e.g. 'RTML-Athena++', 'TIGRESS-ISM', 'FIREbox', 'GIBLE')")
    ap.add_argument("--id", type=str, default=None,
                    help="Object/cloud/box identifier badge (e.g. 'RTML_50pc', 'ISM_1kpc', 'Cloud_12')")
    ap.add_argument("--workers", type=int, default=4,
                    help="Worker processes for WebP lossless atlas encoding and kernel scatter (default 4)")
    args = ap.parse_args(argv)

    if args.input is None and args.preset is None:
        ap.error("Provide either --input <file> or --preset <mixing_layer|ism_kpc>")

    if args.half_pc is None and args.half_kpc is None:
        args.half_pc = 500.0 if args.preset == "ism_kpc" else 25.0

    levels_half_kpc = parse_levels(args.levels, args.half_pc, args.half_kpc)
    ingest_dataset(
        input_path=args.input,
        preset=args.preset,
        out_dir=args.out,
        levels_half_kpc=levels_half_kpc,
        n=args.n,
        nrays=args.nrays,
        sim_name=args.sim,
        subhalo_id=args.id,
        workers=args.workers,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
