# Seeing the Interstellar Medium & Multiphase Mixing Layers (`ism.tanmaysingh.uk`)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Schema: tng50-viz/2](https://img.shields.io/badge/3D_Schema-tng50--viz%2F2-7db2f0.svg)](BYO_SIMULATION.md)
[![Spectra: TNGS v2](https://img.shields.io/badge/Spectra-12--ion_Voigt_TNGS_v2-f39a5c.svg)](BYO_SIMULATION.md)

An interactive, client-side **3D WebGL2 radiative-transfer volume explorer** and **12-ion UV/EUV/X-ray Voigt absorption spectrum laboratory** for the multiphase **Interstellar Medium (ISM)** and **Radiative Turbulent Mixing Layers (RTMLs)** across sub-parsec ($\mathrm{AU}$ to $\mathrm{pc}$) and kiloparsec ($\mathrm{kpc}$) scales.

* **Live Portal:** [https://ism.tanmaysingh.uk](https://ism.tanmaysingh.uk)
* **Companion Galactic Halo / CGM Explorer (`TNG50` & `FIRE-2`):** [https://galaxies.tanmaysingh.uk](https://galaxies.tanmaysingh.uk) ([`tanmay-starslab/cgm-absorber-requirements`](https://github.com/tanmay-starslab/cgm-absorber-requirements))
* **Bring-Your-Own-Simulation CLI Guide:** [`BYO_SIMULATION.md`](BYO_SIMULATION.md)

---

## Key Capabilities

### 1. Sub-Parsec to Kiloparsec 3D WebGL2 Ray-Marcher (`AU` $\to$ `pc` $\to$ `kpc`)
* **Multi-scale nested octree / block hierarchy (`L0..L4`):** Ray-marches up to 5 nested $256^3$ 3D volumes simultaneously in a single GLSL3 fragment pass with exact linear-space 8-corner trilinear interpolation (`texelFetch`), progressive R2 sub-pixel anti-aliasing, and adaptive 60 fps frame-rate governance.
* **Exact physical column integration at any scale:** Whether the box is a $1\ \mathrm{pc}$ turbulent mixing layer ($2000\ \mathrm{AU}$ voxels), a $50\ \mathrm{pc}$ cloud-wind shear interface, or a $1\ \mathrm{kpc}$ stratified ISM fountain patch, integral mode integrates $N_i = \int n_i(\mathbf{x})\,ds$ in physical $\mathrm{cm^{-2}}$ and automatically formats scale bars and publication axes in $\mathrm{AU}$, $\mathrm{pc}$, or $\mathrm{kpc}$.
* **18 physical 3D fields per level:**
  * **Gas state:** Hydrogen number density $n_{\mathrm{H}}\ [\mathrm{cm^{-3}}]$, temperature $T\ [\mathrm{K}]$, metallicity $Z\ [Z_\odot]$
  * **12-ion ionization ladder ($13.6\ \mathrm{eV} \to 739\ \mathrm{eV}$):** $\mathrm{H\,I\ 1216}$, $\mathrm{Mg\,II\ 2796}$, $\mathrm{Si\,II\ 1260}$, $\mathrm{C\,II\ 1335}$, $\mathrm{Si\,III\ 1206}$, $\mathrm{C\,III\ 977}$, $\mathrm{Si\,IV\ 1403}$, $\mathrm{C\,IV\ 1548}$, $\mathrm{N\,V\ 1239}$, $\mathrm{O\,VI\ 1032}$, $\mathrm{Ne\,VIII\ 770}$, $\mathrm{O\,VII\ 21.6\ \mathring{A}}$
  * **5 real-time 3D kinematic views** derived per sample from the 3D velocity vector $(v_x, v_y, v_z)$: line-of-sight velocity $v_{\mathrm{los}}$, line-of-sight velocity dispersion $\sigma_{\mathrm{los}}$, radial velocity $v_r$, azimuthal/shear velocity $v_\phi$, and 3D speed $|\mathbf{v}|$, weighted by total gas ($n_{\mathrm{H}}$) or any selected ion species.
* **1D & 2D Gradient Transfer-Function Volume Rendering:** Isolate intermediate-temperature radiative mixing interfaces ($T \sim 10^{4.8}\text{--}10^{5.5}\ \mathrm{K}$, where $\mathrm{O\,VI}$, $\mathrm{N\,V}$, $\mathrm{C\,IV}$, and $\mathrm{Si\,IV}$ peak) using value $\times$ $|\nabla \log X|$ 2D transfer functions, emissivity weighting ($\propto n_{\mathrm{H}}^\gamma$), and two-sided Blinn–Phong headlight shading.

### 2. Real-Time 12-Ion Voigt Absorption Spectrum Engine
* Click any 3D sightline piercing the mixing layer or ISM volume to inspect its 12-ion absorption spectra degraded live in the browser to any resolving power $R = \lambda/\Delta\lambda \in [5{,}000, 150{,}000]$ and signal-to-noise ratio $\mathrm{S/N} \in [5, 200]$ per resolution element.
* **In-browser Levenberg–Marquardt Voigt fitter (`js/fitworker.js`):** Fits multi-component Voigt profiles $(v_c, \log_{10} N, b)$ through the instrumental line-spread function with $\Delta\mathrm{BIC}$ model selection and saturation flags.
* **Vector Publication Export (`js/export3d.js`):** One-click export of 400 dpi PNG and vector PDF figures with embedded Computer Modern LaTeX math typography and physical $\mathrm{AU} / \mathrm{pc} / \mathrm{kpc}$ axes.

### 3. Bring Your Own Simulation (`python -m tngsight.ingest`)
Convert any private or public `yt`-loadable simulation snapshot/cutout (FIREbox, GIBLE, GLOW, TIGRESS, SILCC, STARFORGE, Athena++, Enzo, RAMSES, AREPO, GIZMO) or any 3D HDF5/NumPy `.npz` cube into a ready-to-serve static web pack in one command:

```bash
python -m tngsight.ingest \
  --input /path/to/simulation_or_grid.hdf5 \
  --out docs/ism/data \
  --half-pc 25 \
  --levels L0,L1,L2 \
  --nrays 120
```

See [`BYO_SIMULATION.md`](BYO_SIMULATION.md) for complete documentation and examples.

---

## Showcase Datasets Included

1. **`RTML_50pc` (`RTML-Athena++` preset, default view):**
   * $50\ \mathrm{pc}$ 3D Radiative Turbulent Mixing Layer ($L_0 = \pm 25\ \mathrm{pc}$, $L_1 = \pm 6.25\ \mathrm{pc}$, $L_2 = \pm 1.5625\ \mathrm{pc}$ down to $0.012\ \mathrm{pc} \approx 2517\ \mathrm{AU}$ resolution) between a cold ISM stream ($T = 10^4\ \mathrm{K}$, $n_{\mathrm{H}} = 1\ \mathrm{cm^{-3}}$, $v_x = -20\ \mathrm{km\,s^{-1}}$) and a hot coronal wind ($T = 10^6\ \mathrm{K}$, $n_{\mathrm{H}} = 10^{-2}\ \mathrm{cm^{-3}}$, $v_x = +140\ \mathrm{km\,s^{-1}}$), resolving fractal Kelvin–Helmholtz rollers, entrained cold droplets, radiative cooling inflow ($v_z \approx -15\text{ to }-35\ \mathrm{km\,s^{-1}}$), and $\mathrm{O\,VI} / \mathrm{C\,IV} / \mathrm{Si\,IV}$ interface layers.
2. **`ISM_1kpc` (`TIGRESS-ISM` preset, selectable in the top picker):**
   * $1\ \mathrm{kpc}$ 3D Stratified Multiphase ISM & Superbubble Fountain Patch ($L_0 = \pm 500\ \mathrm{pc}$, $L_1 = \pm 125\ \mathrm{pc}$, $L_2 = \pm 31.25\ \mathrm{pc}$), resolving a turbulent CNM/WNM midplane, expanding supernova-driven hot superbubble chimneys ($T \sim 10^{6.1}\ \mathrm{K}$, $v_{\mathrm{out}} \sim 130\text{--}145\ \mathrm{km\,s^{-1}}$), and intermediate-ion mixing shells.

---

## License

MIT License — see [`LICENSE`](LICENSE) and [`THIRD_PARTY_LICENSES`](THIRD_PARTY_LICENSES).
