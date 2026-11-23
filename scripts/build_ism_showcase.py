#!/usr/bin/env python3
"""Build the default 3D Radiative Turbulent Mixing Layer (`RTML_50pc`) and Stratified 1 kpc Multiphase ISM
(`ISM_1kpc`) showcase packs under `docs/ism/data/` and render `docs/ism/index.html` + `docs/index.html`."""
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


def render_ism_index(ism_dir: Path, n_rays: int = 120) -> None:
    """Render `docs/ism/index.html` from `docs/ism/_template.html` with the ISM & Mixing-Layer config."""
    tpl = (ism_dir / "_template.html").read_text()
    if MARK not in tpl:
        raise RuntimeError(f"Placeholder {MARK!r} not found in {ism_dir / '_template.html'}")

    cat = json.loads((ism_dir / "data" / "catalog.json").read_text())
    gal_meta = cat["galaxy"]

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
            "base": "data/viz/",
            "fallback_base": None,
            "dataset": "tanmay-starslab/ism-absorber-explorer",
        },
        "js_base": "js/",
        "three": "0.169.0",
        "chat": {"webllm": "0.2.85", "model": "Qwen2.5-1.5B-Instruct-q4f16_1-MLC"},
        "custom": {"r_max_kpc": 0.5, "l_min_kpc": 0.001, "l_max_kpc": 1.0},
        "galaxies": {
            "default": "RTML_50pc",
            "base": "data/gal/",
            "viz_base": None,
            "list": [
                {
                    "id": "RTML_50pc",
                    "sim": "RTML-Athena++",
                    "label": "RTML 50 pc (Mixing Layer)",
                    "Mstar": 1.0e5,
                    "SFR": 0.15,
                    "M200c": 1.0e7,
                    "R200c": 0.025,
                    "central": True,
                    "mw": None,
                    "m31": None,
                    "n": len(cat["entries"]),
                    "viz": True,
                    "src": "trident",
                },
                {
                    "id": "ISM_1kpc",
                    "sim": "TIGRESS-ISM",
                    "label": "ISM 1 kpc (Multiphase Fountain)",
                    "Mstar": 8.0e8,
                    "SFR": 1.25,
                    "M200c": 8.0e10,
                    "R200c": 0.5,
                    "central": True,
                    "mw": True,
                    "m31": None,
                    "n": n_rays,
                    "viz": True,
                    "src": "trident",
                },
            ],
        },
    }

    blob = json.dumps(cfg, allow_nan=False, separators=(",", ":"))
    html = tpl.replace(MARK, blob, 1)
    (ism_dir / "index.html").write_text(html)

    # Also write docs/index.html redirecting cleanly to ism/ so both / and /ism/ work
    root_index = ism_dir.parent / "index.html"
    root_index.write_text(
        '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta http-equiv="refresh" content="0; url=ism/">\n'
        '<title>Seeing the Interstellar Medium &amp; Multiphase Mixing Layers</title>\n'
        '<script>location.replace("ism/" + location.search + location.hash);</script>\n'
        '</head>\n<body><p>Redirecting to <a href="ism/">ISM &amp; Multiphase Mixing Layer Explorer</a>…</p></body>\n</html>\n'
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build the ISM & Turbulent Mixing Layer showcase datasets.")
    ap.add_argument("--levels", type=str, default="L0,L1", help="Levels to build (default L0,L1)")
    ap.add_argument("--n", type=int, default=256, help="Cube resolution per axis (default 256)")
    ap.add_argument("--nrays", type=int, default=120, help="Sightlines per dataset (default 120)")
    ap.add_argument("--workers", type=int, default=8, help="Worker processes (default 8)")
    args = ap.parse_args(argv)

    ism_dir = ROOT / "docs" / "ism"
    # 1. Default dataset: 50 pc Radiative Turbulent Mixing Layer (L0 = ±25 pc, L1 = ±6.25 pc, L2 = ±1.5625 pc)
    lv_rtml = parse_levels(args.levels, half_pc=25.0, half_kpc=None)
    ingest_dataset(
        preset="mixing_layer",
        out_dir=ism_dir / "data",
        levels_half_kpc=lv_rtml,
        n=args.n,
        nrays=args.nrays,
        sim_name="RTML-Athena++",
        subhalo_id="RTML_50pc",
        workers=args.workers,
    )
    VP.validate_manifest_v2(
        json.loads((ism_dir / "data" / "viz" / "manifest_v2.json").read_text()),
        base_dir=ism_dir / "data" / "viz",
    )

    # 2. Companion dataset: 1 kpc Stratified Multiphase ISM & Superbubble Fountain Patch (L0 = ±500 pc, L1 = ±125 pc)
    lv_ism = parse_levels(args.levels, half_pc=500.0, half_kpc=None)
    gal_ism_dir = ism_dir / "data" / "gal" / "ISM_1kpc"
    ingest_dataset(
        preset="ism_kpc",
        out_dir=gal_ism_dir,
        levels_half_kpc=lv_ism,
        n=args.n,
        nrays=args.nrays,
        sim_name="TIGRESS-ISM",
        subhalo_id="ISM_1kpc",
        workers=args.workers,
    )
    VP.validate_manifest_v2(
        json.loads((gal_ism_dir / "viz" / "manifest_v2.json").read_text()),
        base_dir=gal_ism_dir / "viz",
    )

    render_ism_index(ism_dir, n_rays=args.nrays)
    print(f"[build_ism_showcase] Successfully built & validated RTML_50pc + ISM_1kpc and rendered {ism_dir / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
