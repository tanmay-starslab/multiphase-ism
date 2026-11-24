#!/usr/bin/env python3
"""Upload all 6 flagship 5-level (L0..L4 gas + S0..S4 starlight/dust) 3D volume directories
and sightline packs from `build_5lvl/` to Hugging Face dataset `tanmay-starslab/ism-atlas-3d`,
and sync lightweight web sightline packs + manifests into `docs/ism/data/`."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

from huggingface_hub import HfApi

ROOT = Path(__file__).resolve().parents[1]
REPO_ID = "tanmay-starslab/ism-atlas-3d"
VOLUMES = [
    "RTML_50pc",
    "CloudCrush_200pc",
    "GMC_STARFORGE_20pc",
    "Superbubble_500pc",
    "ShatteringMist_10pc",
    "ISM_1kpc",
]

HF_README = """---
license: mit
task_categories:
  - other
tags:
  - astrophysics
  - interstellar-medium
  - turbulent-mixing-layers
  - cloud-crushing
  - starforge
  - superbubble
  - thermal-instability
  - spectroscopy
  - webgl2
  - 3d-volume-rendering
pretty_name: Multiphase ISM & Sub-ISM 3D Octree & 12-Ion Voigt Sightline Atlas (31 AU to 1 kpc)
---

# Multiphase ISM & Sub-ISM 3D Octree & 12-Ion Voigt Sightline Atlas (`tanmay-starslab/ism-atlas-3d`)

Public 5-level ($256^3 \\times 5$ nested octree levels: `L0..L4` gas/ions/velocities + `S0..S4` starlight/dust) 3D WebGL2 volume dataset and 720-sightline 12-ion FP64 Voigt absorption library backing **[ism.tanmaysingh.uk](https://ism.tanmaysingh.uk)** ([tanmay-starslab/multiphase-ism](https://github.com/tanmay-starslab/multiphase-ism)).

## 6 Flagship Modern ISM & Sub-ISM Regimes ($31\\text{ AU} \\rightarrow 1\\text{ kpc}$, 720 Sightlines)

| Volume ID | `simKey` | Astrophysical Regime | Box Extent ($2H_0$) | 5-Level Octree (`L0` $\\rightarrow$ `L4`, `S0` $\\rightarrow$ `S4`) | Finest Voxel (`L4`) | Sightlines |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: |
| **`RTML_50pc`** | `rtml` | Radiative Turbulent Mixing Layer (Athena++ / Fielding+2020 / Tan+2021) | $50\\text{ pc}$ | $\\pm 25\\text{ pc} \\rightarrow \\pm 0.0977\\text{ pc}$ | $0.00076\\text{ pc}$ ($157\\text{ AU}$) | 120 |
| **`CloudCrush_200pc`** | `cloudcrush` | Radiative Cloud Crushing & Cometary Tail (Gronke-Oh / Cholla) | $200\\text{ pc}$ | $\\pm 100\\text{ pc} \\rightarrow \\pm 0.3906\\text{ pc}$ | $0.00305\\text{ pc}$ ($630\\text{ AU}$) | 120 |
| **`GMC_STARFORGE_20pc`** | `starforge` | Magnetized GMC & Blister H II Bubble (STARFORGE / SILCC) | $20\\text{ pc}$ | $\\pm 10\\text{ pc} \\rightarrow \\pm 0.0391\\text{ pc}$ | $0.000305\\text{ pc}$ ($62\\text{ AU}$) | 120 |
| **`Superbubble_500pc`** | `superbubble` | Clustered-SN Superbubble Chimney Breakout (TIGRESS / PRFM) | $500\\text{ pc}$ | $\\pm 250\\text{ pc} \\rightarrow \\pm 0.9766\\text{ pc}$ | $0.00763\\text{ pc}$ ($1570\\text{ AU}$) | 120 |
| **`ShatteringMist_10pc`** | `shattering` | Isobaric Thermal Instability Shattering Mist (McCourt-Gronke $c_s t_{\\text{cool}}$) | $10\\text{ pc}$ | $\\pm 5\\text{ pc} \\rightarrow \\pm 0.0195\\text{ pc}$ | $0.000153\\text{ pc}$ ($31\\text{ AU}$) | 120 |
| **`ISM_1kpc`** | `ism` | Stratified Multiphase ISM & Superbubble Fountain (TIGRESS-ISM) | $1\\text{ kpc}$ | $\\pm 500\\text{ pc} \\rightarrow \\pm 1.9531\\text{ pc}$ | $0.01526\\text{ pc}$ ($3148\\text{ AU}$) | 120 |

## Directory Structure per Volume (`<Volume_ID>/`)
- `viz/manifest_v2.json` — `tng50-viz/2` 5-level gas (`L0..L4`) + 5-level starlight/dust (`S0..S4`) manifest with SHA-256 checksums, log quantization ranges, and parsec/AU metadata
- `viz/gas/L{0..4}/` — 18 fields per level (`nH`, `T`, `Z`, `HI`, `MgII`, `SiII`, `CII`, `SiIII`, `CIII`, `SiIV`, `CIV`, `NV`, `OVI`, `NeVIII`, `OVII`, `vx`, `vy`, `vz`) stored as 4 lossless $1024\\times 4096$ WebP z-slice atlases (`.p0..3.webp`) + raw `.u8.gz`
- `viz/stars/S{0..4}/` — 4-channel RGBA starlight ($i, r, g$) + V-band dust extinction ($\\kappa_V$) cubes (`.p0..3.webp` + `.rgba.gz`)
- `viz/luts/luts.json` — colormap lookup tables
- `catalog.json`, `components.json`, `galaxy.json`, `s/*.bin` — 120-ray 12-ion `TNGS` v2 sightline pack
"""


def get_token() -> str:
    if os.environ.get("HF_TOKEN"):
        return os.environ["HF_TOKEN"].strip()
    tok_path = Path("/scratch/tsingh65/cgm_work/cache/hf/token")
    if tok_path.is_file():
        return tok_path.read_text().strip()
    raise RuntimeError("HF_TOKEN not found")


def sync_web_packs_to_docs(build_root: Path, docs_ism_dir: Path, volumes: list[str] | None = None) -> None:
    """Copy lightweight web packs (catalog.json, components.json, galaxy.json, s/*.bin, viz/manifest_v2.json, viz/luts/luts.json)
    into `docs/ism/data/gal/<id>/` (and `docs/ism/data/` for default `RTML_50pc`), removing heavy local cube directories
    from git so all 5 levels stream cleanly from Hugging Face `tanmay-starslab/ism-atlas-3d`."""
    data_dir = docs_ism_dir / "data"
    for vid in (volumes or VOLUMES):
        src = build_root / vid
        if not (src / "catalog.json").is_file() or not (src / "viz" / "manifest_v2.json").is_file():
            raise FileNotFoundError(f"Missing built volume in {src}")

        targets = [data_dir / "gal" / vid]
        if vid == "RTML_50pc":
            targets.append(data_dir)

        for dst in targets:
            dst.mkdir(parents=True, exist_ok=True)
            for fn in ("catalog.json", "components.json", "galaxy.json"):
                shutil.copy2(src / fn, dst / fn)
            dst_s = dst / "s"
            if dst_s.exists():
                shutil.rmtree(dst_s)
            shutil.copytree(src / "s", dst_s)

            dst_viz = dst / "viz"
            dst_viz.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src / "viz" / "manifest_v2.json", dst_viz / "manifest_v2.json")
            if (src / "viz" / "luts" / "luts.json").is_file():
                (dst_viz / "luts").mkdir(parents=True, exist_ok=True)
                shutil.copy2(src / "viz" / "luts" / "luts.json", dst_viz / "luts" / "luts.json")

            # Remove heavy local gas/stars cubes from docs/ism/data so git stays lightweight
            for heavy in ("gas", "stars"):
                if (dst_viz / heavy).exists():
                    shutil.rmtree(dst_viz / heavy)
        print(f"[sync] Synced web pack + 5-level manifest for {vid} -> {', '.join(str(t) for t in targets)}", flush=True)


def upload_to_hf(build_root: Path, volumes: list[str] | None = None) -> None:
    token = get_token()
    api = HfApi(token=token)
    api.create_repo(repo_id=REPO_ID, repo_type="dataset", private=False, exist_ok=True)

    readme_path = build_root / "README.md"
    readme_path.write_text(HF_README)
    api.upload_file(
        path_or_fileobj=str(readme_path),
        path_in_repo="README.md",
        repo_id=REPO_ID,
        repo_type="dataset",
        commit_message="Update dataset card for 6 flagship ISM & sub-ISM 5-level 3D volumes",
    )

    vols = volumes or VOLUMES
    for vid in vols:
        vdir = build_root / vid
        if not (vdir / "viz" / "manifest_v2.json").is_file():
            raise FileNotFoundError(f"Missing {vdir / 'viz' / 'manifest_v2.json'}")
        print(f"[hf] Uploading {vid} (5 levels L0..L4 + S0..S4 + 120 sightlines) to {REPO_ID}/{vid}...", flush=True)
        api.upload_folder(
            repo_id=REPO_ID,
            repo_type="dataset",
            folder_path=str(vdir),
            path_in_repo=vid,
            commit_message=f"Upload 5-level (L0..L4 + S0..S4) 3D volume and 120-ray library for {vid}",
        )
        print(f"[hf] Finished uploading {vid} -> https://huggingface.co/datasets/{REPO_ID}/tree/main/{vid}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-root", type=Path, default=ROOT / "build_5lvl")
    ap.add_argument("--volumes", type=str, default=None, help="Comma-separated subset of volumes to upload/sync")
    ap.add_argument("--sync-only", action="store_true")
    ap.add_argument("--upload-only", action="store_true")
    args = ap.parse_args()

    vols = [x.strip() for x in args.volumes.split(",") if x.strip()] if args.volumes else None
    if not args.upload_only:
        sync_web_packs_to_docs(args.build_root, ROOT / "docs" / "ism", volumes=vols)
    if not args.sync_only:
        upload_to_hf(args.build_root, volumes=vols)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
