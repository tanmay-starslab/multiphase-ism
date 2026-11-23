# Bring-Your-Own-Simulation (`python -m tngsight.ingest`) Guide

`multiphase-ism` includes a universal one-command ingestion pipeline (`tngsight.ingest`) that converts any 3D astrophysical hydro/MHD simulation snapshot, cutout, or uniform grid into a self-contained static 3D WebGL2 + 12-ion Voigt absorption sightline web pack (`tng50-viz/2` + `TNGS` v2).

Because the entire explorer runs client-side in the browser (WebGL2 3D radiative transfer + 60 fps 1D/2D transfer functions + Web Worker Levenberg–Marquardt Voigt fitter), researchers working with **private or non-public simulations** (such as **FIREbox**, **GIBLE**, **GLOW**, **TIGRESS**, **SILCC**, **STARFORGE**, or custom **Athena++ / Enzo / RAMSES / AREPO / GIZMO** runs) can visualize and fit their own volumes locally on their laptop or workstation without uploading proprietary simulation snapshots to any external server.

---

## 1. Quick Start (One Command)

### A. Convert any `yt`-loadable simulation snapshot or cutout
Supports **Athena++** (`.athdf`), **AREPO / IllustrisTNG / GIBLE / GLOW** (`.hdf5`), **GIZMO / Gadget / FIRE-2 / FIREbox / STARFORGE** (`.hdf5`), **Enzo**, **RAMSES**, **FLASH / SILCC**, and **TIGRESS** VTK/Athena++ outputs:

```bash
python -m tngsight.ingest \
  --input /path/to/snapshot_or_cutout.hdf5 \
  --out docs/ism/data \
  --half-pc 25 \
  --levels L0,L1,L2 \
  --n 256 \
  --nrays 120 \
  --sim "Athena++-RTML" \
  --id "KHI_Cloud_50pc"
```

Then serve `docs/ism/` locally with any static HTTP server:

```bash
python -m http.server 8080 --directory docs/ism
# Open http://localhost:8080/ in any WebGL2 browser
```

### B. Convert a raw 3D NumPy (`.npz`) or HDF5 (`.h5`) grid (Zero `yt` dependency)
If you have a 3D Cartesian grid or particle table in Python, save it to `.npz` or `.h5` with fields:
- `nH`: 3D float array `(N, N, N)` of hydrogen number density in $\mathrm{cm^{-3}}$
- `T`: 3D float array `(N, N, N)` of gas temperature in $\mathrm{K}$
- `Z` *(optional)*: 3D float array `(N, N, N)` of metallicity in $Z_\odot$ (defaults to $1.0\,Z_\odot$)
- `vx`, `vy`, `vz` *(optional)*: 3D float arrays `(N, N, N)` of velocity components in $\mathrm{km\,s^{-1}}$
- Any precomputed ion number density arrays (`HI`, `MgII`, `SiII`, `CII`, `SiIII`, `CIII`, `SiIV`, `CIV`, `NV`, `OVI`, `NeVIII`, `OVII` in $\mathrm{cm^{-3}}$); any omitted ions are computed automatically from $(n_{\mathrm{H}}, T, Z)$ using CIE (Gnat & Sternberg 2007) + Haardt & Madau (2012) UV photoionization + Rahmati et al. (2013) self-shielding.

```python
import numpy as np

# Example: save a 256^3 custom MHD mixing-layer cube (50 pc box -> half_pc = 25)
np.savez_compressed(
    "my_mixing_layer.npz",
    nH=nH_cube,   # [cm^-3]
    T=T_cube,     # [K]
    Z=Z_cube,     # [Z_sun]
    vx=vx_cube,   # [km/s]
    vy=vy_cube,   # [km/s]
    vz=vz_cube,   # [km/s]
)
```

Ingest into a full 3D WebGL2 + 120-ray Voigt spectrum pack:

```bash
python -m tngsight.ingest \
  --input my_mixing_layer.npz \
  --out docs/ism/data \
  --half-pc 25 \
  --levels L0,L1 \
  --sim "Custom-MHD" \
  --id "Layer_25pc"
```

---

## 2. Multi-Scale Coordinate & Level Hierarchy (`AU` $\to$ `pc` $\to$ `kpc`)

`manifest_v2.json` (`tng50-viz/2` schema) stores physical coordinates in canonical kiloparsecs (`half_kpc`, `voxel_kpc`) together with parsec metadata (`half_pc`, `voxel_pc`), and the WebGL2 viewer (`viewer3d.js`, `raymarch.js`, `export3d.js`) automatically adapts its camera frustum, ray-march step $ds$, volume optical depth $\kappa$, and physical scale bar across **9 orders of magnitude**:

| Physical Regime | Example `--half-pc` | `L0` Box Width | `L2` Voxel Size ($256^3$) | Scale Bar Units |
|---|---|---|---|---|
| Sub-parsec cloudlet / shock filament | `--half-pc 0.5` | $1\ \mathrm{pc}$ | $4.9\times 10^{-4}\ \mathrm{pc} \approx 101\ \mathrm{AU}$ | `AU` / `pc` |
| Radiative Turbulent Mixing Layer (RTML) | `--half-pc 25` | $50\ \mathrm{pc}$ | $0.0122\ \mathrm{pc} \approx 2517\ \mathrm{AU}$ | `pc` / `AU` |
| Giant Molecular Cloud / STARFORGE core | `--half-pc 100` | $200\ \mathrm{pc}$ | $0.0488\ \mathrm{pc}$ | `pc` |
| Stratified Multiphase ISM Patch (TIGRESS / SILCC) | `--half-pc 500` | $1\ \mathrm{kpc}$ | $0.244\ \mathrm{pc}$ | `pc` / `kpc` |
| Galactic Halo / CGM (TNG50 / FIRE-2 / GIBLE) | `--half-kpc 400` | $800\ \mathrm{kpc}$ | $195\ \mathrm{pc}$ ($12\ \mathrm{pc}$ at `L4`) | `kpc` / `pc` |

Every exact integral ray march step multiplies the local density $n_i(\mathbf{x})\ [\mathrm{cm^{-3}}]$ by the true physical path element $ds \cdot (3.08567758\times 10^{21}\ \mathrm{cm\,kpc^{-1}})$, so projected column densities $\log_{10} N_i\ [\mathrm{cm^{-2}}]$ are physically exact at every scale from $100\ \mathrm{AU}$ to $1\ \mathrm{Mpc}$.

---

## 3. Built-in 3D Physical Presets

To test or benchmark the explorer without external snapshot files, `tngsight.ingest` provides two physically realistic 3D multiphase presets:

1. **Radiative Turbulent Mixing Layer (`--preset mixing_layer`)**:
   - Models a cold ISM stream ($T = 10^4\ \mathrm{K}$, $n_{\mathrm{H}} = 1.0\ \mathrm{cm^{-3}}$, $Z = 1.0\,Z_\odot$, $v_x = -20\ \mathrm{km\,s^{-1}}$) shearing against a hot galactic wind/coronal medium ($T = 10^6\ \mathrm{K}$, $n_{\mathrm{H}} = 10^{-2}\ \mathrm{cm^{-3}}$, $Z = 0.35\,Z_\odot$, $v_x = +140\ \mathrm{km\,s^{-1}}$, $\Delta v = 160\ \mathrm{km\,s^{-1}}$) across a fractal Kelvin–Helmholtz interface with entrained cold droplets and a $T \sim 10^{4.9}\text{--}10^{5.5}\ \mathrm{K}$ radiative cooling boundary layer (rich in $\mathrm{O\,VI}$, $\mathrm{N\,V}$, $\mathrm{C\,IV}$, $\mathrm{Si\,IV}$, $\mathrm{C\,III}$) and radiative cooling inflow $v_z \approx -15\text{ to }-35\ \mathrm{km\,s^{-1}}$ (cf. Ji, Oh & Masterson 2019; Fielding et al. 2020; Tan, Oh & Gronke 2021).
   ```bash
   python -m tngsight.ingest --preset mixing_layer --out docs/ism/data --half-pc 25 --levels L0,L1,L2
   ```

2. **Stratified 1 kpc Multiphase ISM & Superbubble Fountain Patch (`--preset ism_kpc`)**:
   - Models a $1\ \mathrm{kpc}$ stratified turbulent ISM patch ($L_0 = \pm 500\ \mathrm{pc}$, $L_1 = \pm 125\ \mathrm{pc}$, $L_2 = \pm 31.25\ \mathrm{pc}$) with a cold neutral/molecular midplane, warm neutral/ionized filaments, and expanding supernova-driven hot superbubble chimneys ($T \sim 10^{6.1}\ \mathrm{K}$, $v_{\mathrm{out}} \sim 130\text{--}145\ \mathrm{km\,s^{-1}}$) coated in intermediate-ion mixing shells (cf. Kim & Ostriker 2017 TIGRESS; Walch et al. 2015 SILCC).
   ```bash
   python -m tngsight.ingest --preset ism_kpc --out docs/ism/data/gal/ISM_1kpc --half-pc 500 --levels L0,L1,L2
   ```

---

## 4. Output Directory Contract

Running `python -m tngsight.ingest --out <DIR>` produces:

```text
<DIR>/
├── catalog.json              # Sightline catalogue v2 (endpoints p0/p1, logN, AOD columns, scores, ranks)
├── components.json           # Fitted Voigt component catalogue (v, dv, logN, dlogN, b, db, EW, flag)
├── galaxy.json               # Simulation & coordinate frame metadata
├── s/
│   ├── face_i000.bin         # TNGS v2 binary spectra (12 transitions per sightline)
│   ├── edge_i000.bin
│   └── rand_s0001.bin
└── viz/
    ├── manifest_v2.json      # Validated tng50-viz/2 3D volume manifest (SHA-256 + byte sizes)
    ├── luts/luts.json        # Perceptually uniform colour maps
    └── gas/
        ├── manifest_gas.json
        ├── L0/               # 18 fields × (4 lossless WebP z-slice atlases + 1 .u8.gz cube)
        ├── L1/
        └── L2/
```

Every generated `manifest_v2.json` is automatically verified against all cube files on disk via `python -m tngsight.vizpack validate <DIR>/viz`.
