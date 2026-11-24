#!/usr/bin/env python3
"""Build all 6 flagship modern ISM & sub-ISM 5-level (L0..L4 gas + S0..S4 starlight/dust) simulation packs
and render `docs/ism/index.html` + `docs/ism/data/galaxies.json` + `docs/index.html` for `ism.tanmaysingh.uk`:
  1. `RTML_50pc` (`simKey: "rtml"`, 50 pc -> 157 AU Radiative Turbulent Mixing Layer, Athena++ style)
  2. `CloudCrush_200pc` (`simKey: "cloudcrush"`, 200 pc -> 630 AU Radiative Cloud Crushing & Cometary Tail, Gronke-Oh / Cholla)
  3. `GMC_STARFORGE_20pc` (`simKey: "starforge"`, 20 pc -> 62 AU Magnetized GMC & Blister H II Bubble, STARFORGE / SILCC)
  4. `Superbubble_500pc` (`simKey: "superbubble"`, 500 pc -> 1570 AU Clustered-SN X-ray Superbubble Breakout, TIGRESS / PRFM)
  5. `ShatteringMist_10pc` (`simKey: "shattering"`, 10 pc -> 31 AU Isobaric Thermal Instability Shattering Mist, McCourt-Gronke)
  6. `ISM_1kpc` (`simKey: "ism"`, 1 kpc -> 3148 AU Stratified Multiphase ISM & Superbubble Fountain, TIGRESS style)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tngsight import ions as I
from tngsight import vizpack as VP
from tngsight.ingest import ingest_dataset, parse_levels

MARK = "/*__DATA__*/null"
HF_DATASET = "tanmay-starslab/ism-atlas-3d"
HF_VIZ_BASE = f"https://huggingface.co/datasets/{HF_DATASET}/resolve/main/{{id}}/viz/"

VOLUMES_SPEC = [
    {
        "id": "RTML_50pc",
        "preset": "mixing_layer",
        "sim": "RTML-Athena++",
        "simKey": "rtml",
        "label": "50 pc RTML",
        "half_pc": 25.0,
        "Mstar": 1.0e5,
        "SFR": 0.15,
        "M200c": 1.0e7,
        "R200c": 0.025,
        "central": True,
        "mw": None,
        "m31": None,
        "n": 120,
        "viz": True,
        "src": "trident",
    },
    {
        "id": "CloudCrush_200pc",
        "preset": "cloud_crush",
        "sim": "Gronke-Oh / Cholla Cloud Crushing",
        "simKey": "cloudcrush",
        "label": "200 pc Cloud Crush",
        "half_pc": 100.0,
        "Mstar": 3.5e5,
        "SFR": 0.42,
        "M200c": 4.8e7,
        "R200c": 0.100,
        "central": True,
        "mw": None,
        "m31": None,
        "n": 120,
        "viz": True,
        "src": "trident",
    },
    {
        "id": "GMC_STARFORGE_20pc",
        "preset": "starforge_gmc",
        "sim": "STARFORGE / SILCC GMC & H II Bubble",
        "simKey": "starforge",
        "label": "20 pc STARFORGE GMC",
        "half_pc": 10.0,
        "Mstar": 2.4e4,
        "SFR": 0.085,
        "M200c": 2.0e5,
        "R200c": 0.010,
        "central": True,
        "mw": None,
        "m31": None,
        "n": 120,
        "viz": True,
        "src": "trident",
    },
    {
        "id": "Superbubble_500pc",
        "preset": "superbubble",
        "sim": "TIGRESS / PRFM Clustered-SN Superbubble",
        "simKey": "superbubble",
        "label": "500 pc Superbubble",
        "half_pc": 250.0,
        "Mstar": 1.8e7,
        "SFR": 1.25,
        "M200c": 3.2e8,
        "R200c": 0.250,
        "central": True,
        "mw": True,
        "m31": None,
        "n": 120,
        "viz": True,
        "src": "trident",
    },
    {
        "id": "ShatteringMist_10pc",
        "preset": "shattering_mist",
        "sim": "McCourt-Gronke Thermal Instability Mist",
        "simKey": "shattering",
        "label": "10 pc Shattering Mist",
        "half_pc": 5.0,
        "Mstar": 1.2e3,
        "SFR": 0.012,
        "M200c": 4.5e4,
        "R200c": 0.005,
        "central": True,
        "mw": None,
        "m31": None,
        "n": 120,
        "viz": True,
        "src": "trident",
    },
    {
        "id": "ISM_1kpc",
        "preset": "ism_kpc",
        "sim": "TIGRESS-ISM",
        "simKey": "ism",
        "label": "1 kpc ISM",
        "half_pc": 500.0,
        "Mstar": 8.0e8,
        "SFR": 1.25,
        "M200c": 8.0e10,
        "R200c": 0.500,
        "central": True,
        "mw": True,
        "m31": None,
        "n": 120,
        "viz": True,
        "src": "trident",
    },
]


def galaxies_index_payload(n_rays: int = 120) -> dict:
    return {
        "default": "RTML_50pc",
        "base": "data/gal/",
        "viz_base": HF_VIZ_BASE,
        "list": [
            {**{k: v for k, v in spec.items() if k not in ("preset", "half_pc")}, "n": int(n_rays)}
            for spec in VOLUMES_SPEC
        ],
    }


def render_ism_index(ism_dir: Path, n_rays: int = 120) -> None:
    """Render `docs/ism/index.html` from `docs/ism/_template.html` and write `docs/ism/data/galaxies.json`."""
    tpl = (ism_dir / "_template.html").read_text()
    if MARK not in tpl:
        raise RuntimeError(f"Placeholder {MARK!r} not found in {ism_dir / '_template.html'}")

    cat = json.loads((ism_dir / "data" / "catalog.json").read_text())
    gal_meta = cat["galaxy"]
    gix_payload = galaxies_index_payload(n_rays=len(cat["entries"]) or n_rays)
    (ism_dir / "data" / "galaxies.json").write_text(json.dumps(gix_payload, indent=2, allow_nan=False) + "\n")

    cfg = {
        "version": 3,
        "catalog_url": "data/catalog.json",
        "bin_dir": "data/s/",
        "components_url": "data/components.json",
        "n_static": len(cat["entries"]),
        "backends": [],
        "galaxy": gal_meta,
        "ions": I.page_ions(),
        "bands": I.BANDS,
        "snr": I.SNR,
        "defaults": {**I.DEFAULTS, "layout": "1"},
        "grid": I.GRID,
        "registry_tag": I.registry_tag(),
        "fitter": {"worker": "js/fitworker.js", "budget_ms": 300, "cap_ms": 1500, "version": "fit/1"},
        "frame": {
            "Lhat": gal_meta.get("Lhat", [0.0, 0.0, 1.0]),
            "xhat": gal_meta.get("xhat", [1.0, 0.0, 0.0]),
            "yhat": gal_meta.get("yhat", [0.0, 1.0, 0.0]),
            "R200c_kpc": gal_meta.get("R200c_kpc", 0.025),
        },
        "viz": {
            "manifest_url": "data/viz/manifest_v2.json",
            "base": f"https://huggingface.co/datasets/{HF_DATASET}/resolve/main/RTML_50pc/viz/",
            "fallback_base": "data/viz/",
            "dataset": HF_DATASET,
        },
        "js_base": "js/",
        "three": "0.169.0",
        "chat": {"webllm": "0.2.85", "model": "Qwen2.5-1.5B-Instruct-q4f16_1-MLC"},
        "custom": {"r_max_kpc": 0.5, "l_min_kpc": 0.001, "l_max_kpc": 1.0},
        "galaxies": gix_payload,
    }

    blob = json.dumps(cfg, allow_nan=False, separators=(",", ":"))
    html = tpl.replace(MARK, blob, 1)
    (ism_dir / "index.html").write_text(html)

    root_index = ism_dir.parent / "index.html"
    root_index.write_text(
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta http-equiv="refresh" content="0; url=ism/">\n'
        '<title>Seeing the Interstellar Medium &amp; Multiphase Mixing Layers</title>\n'
        '<script>location.replace("ism/" + location.search + location.hash);</script>\n'
        '</head>\n<body><p>Redirecting to <a href="ism/">ISM &amp; Multiphase Mixing Layer Explorer</a>…</p></body>\n</html>\n'
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build the 6 flagship ISM & sub-ISM simulation regimes.")
    ap.add_argument("--render-only", action="store_true", help="Only render docs/ism/index.html and galaxies.json")
    ap.add_argument("--levels", type=str, default="L0,L1,L2,L3,L4", help="Levels to build (default L0,L1,L2,L3,L4)")
    ap.add_argument("--n", type=int, default=256, help="Cube resolution per axis (default 256)")
    ap.add_argument("--nrays", type=int, default=120, help="Sightlines per dataset (default 120)")
    ap.add_argument("--workers", type=int, default=8, help="Worker processes (default 8)")
    ap.add_argument("--only", type=str, default=None, help="Build only one volume id (e.g. CloudCrush_200pc)")
    ap.add_argument("--out-root", type=Path, default=None, help="Optional scratch build root directory")
    args = ap.parse_args(argv)

    ism_dir = ROOT / "docs" / "ism"
    if not args.render_only:
        for spec in VOLUMES_SPEC:
            if args.only and args.only not in (spec["id"], spec["simKey"], spec["preset"]):
                continue
            lv = parse_levels(args.levels, half_pc=spec["half_pc"], half_kpc=None)
            out_dir = (args.out_root / spec["id"]) if args.out_root else (
                (ism_dir / "data") if spec["id"] == "RTML_50pc" else (ism_dir / "data" / "gal" / spec["id"])
            )
            ingest_dataset(
                preset=spec["preset"],
                out_dir=out_dir,
                levels_half_kpc=lv,
                n=args.n,
                nrays=args.nrays,
                sim_name=spec["sim"],
                subhalo_id=spec["id"],
                workers=args.workers,
            )
            VP.validate_manifest_v2(
                json.loads((out_dir / "viz" / "manifest_v2.json").read_text()),
                base_dir=out_dir / "viz",
            )

    if (ism_dir / "data" / "catalog.json").is_file():
        render_ism_index(ism_dir, n_rays=args.nrays)
        print(f"[build_ism_showcase] Rendered {ism_dir / 'index.html'} & {ism_dir / 'data' / 'galaxies.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
