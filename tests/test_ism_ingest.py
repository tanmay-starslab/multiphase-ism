"""Tests for the ISM & Mixing-Layer Bring-Your-Own-Simulation (`tngsight.ingest`) pipeline and parsec/sub-parsec
`manifest_v2.json` validation."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from tngsight import ingest as ING
from tngsight import library as L
from tngsight import vizpack as VP


def test_ion_densities_multiphase_ladder():
    """Verify that H I peaks in cold self-shielded gas (10^4 K), C IV / O VI peak in the T ~ 10^5-10^5.5 K
    radiative mixing regime, and O VII peaks in the T ~ 10^6.2 K hot phase."""
    nH = np.array([1.0, 0.05, 0.01])
    T = np.array([1.0e4, 2.8e5, 1.5e6])
    Z = np.array([1.0, 0.7, 0.5])
    ions = ING.compute_ion_densities(nH, T, Z)
    assert set(ions.keys()) == set(ING.I.KEYS)
    # H I fraction is highest in the cold 10^4 K cell
    assert ions["HI"][0] / nH[0] > 0.5
    assert ions["HI"][2] / nH[2] < 1e-4
    # O VI and C IV peak in the T = 2.8e5 K mixing-layer cell
    assert ions["OVI"][1] / nH[1] > ions["OVI"][0] / nH[0]
    assert ions["OVI"][1] / nH[1] > ions["OVI"][2] / nH[2]
    # O VII peaks in the T = 1.5e6 K hot coronal cell
    assert ions["OVII"][2] / nH[2] > ions["OVII"][0] / nH[0]


def test_byo_simulation_npz_subparsec_ingest(tmp_path: Path):
    """Verify that an external 3D .npz simulation grid at sub-parsec scale (half_pc = 0.5 pc -> 1 pc box)
    ingests into a validated manifest_v2.json + TNGS v2 sightline pack."""
    n = 64
    rng = np.random.default_rng(7)
    nH = 10.0 ** rng.uniform(-2.0, 1.0, size=(n, n, n))
    T = 10.0 ** rng.uniform(3.9, 6.1, size=(n, n, n))
    Z = np.full((n, n, n), 1.0)
    vx = rng.normal(0.0, 35.0, size=(n, n, n))
    vy = rng.normal(0.0, 20.0, size=(n, n, n))
    vz = rng.normal(-15.0, 20.0, size=(n, n, n))

    npz_path = tmp_path / "custom_cloud_1pc.npz"
    np.savez(npz_path, nH=nH, T=T, Z=Z, vx=vx, vy=vy, vz=vz)

    out_dir = tmp_path / "web_pack"
    levels = ING.parse_levels("L0,L1", half_pc=0.5, half_kpc=None)
    assert abs(levels["L0"] - 0.0005) < 1e-12
    assert abs(levels["L1"] - 0.000125) < 1e-12

    res = ING.ingest_dataset(
        input_path=npz_path,
        out_dir=out_dir,
        levels_half_kpc=levels,
        n=64,
        nrays=12,
        sim_name="BYO-Athena++",
        subhalo_id="Cloud_1pc",
        workers=2,
    )
    man_path = Path(res["manifest"])
    assert man_path.is_file()
    man = json.loads(man_path.read_text())
    VP.validate_manifest_v2(man, base_dir=man_path.parent)
    assert man["levels"]["L0"]["half_pc"] == 0.5
    assert man["levels"]["L1"]["half_pc"] == 0.125
    assert man["stars"]["levels"] == []

    cat = json.loads((out_dir / "catalog.json").read_text())
    assert len(cat["entries"]) == 12
    first_id = cat["entries"][0]["id"]
    bin_blob = (out_dir / "s" / f"{first_id}.bin").read_bytes()
    dec = L.decode_web(bin_blob)
    assert dec["version"] == 2
    assert len(dec["lines"]) == 12
