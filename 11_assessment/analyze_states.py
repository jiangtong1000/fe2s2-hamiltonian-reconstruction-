"""Common density diagnostics and historical comparison; no orbital selection."""

from pathlib import Path as _Path
import sys as _sys

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "_shared"))
from workflow_paths import DATA_ROOT, CODE_ROOT
import argparse
import csv
import json
import subprocess
import sys
import numpy as np

ROOT = DATA_ROOT


def analyze(path):
    dm = np.load(path / "rdm1.npy")
    energies = json.loads((path / "energies.json").read_text())
    convergence = json.loads((path / "convergence.json").read_text())
    assert (
        abs(energies["energy_total"] - convergence[-1]["energy_total"]) < 1e-7
    ), "State files and convergence history are from different stages"
    latest = np.load(path / f"rdm1_stage{convergence[-1]['stage']}.npy")
    assert (
        abs(dm - latest).max() < 1e-8
    ), "Final density has not caught up with the latest convergence stage"
    occ, u = np.linalg.eigh((dm + dm.T) / 2)
    occ = occ[::-1]
    u = u[:, ::-1]
    assert abs(np.trace(dm) - energies["nelec"]) < 1e-7
    assert abs(dm - dm.T).max() < 1e-8
    assert occ.min() > -1e-8 and occ.max() < 2 + 1e-8
    change = None
    subspace_change = None
    energy_change = None
    if len(convergence) > 1:
        prev = np.load(path / f"rdm1_stage{convergence[-2]['stage']}.npy")
        old, old_u = np.linalg.eigh((prev + prev.T) / 2)
        change = float(abs(occ - old[::-1]).max())
        singular = np.linalg.svd(old_u[:, -6:].T @ u[:, :6], compute_uv=False)
        subspace_change = float(np.sqrt(max(0.0, 1 - singular.min() ** 2)))
        energy_change = float(
            convergence[-1]["energy_total"] - convergence[-2]["energy_total"]
        )
    result = dict(
        norb=len(occ),
        nelec=energies["nelec"],
        occupations=occ.tolist(),
        density_trace=float(np.trace(dm)),
        hermiticity_max_error=float(abs(dm - dm.T).max()),
        energy_total=energies["energy_total"],
        energy_active=energies["energy_active"],
        ecore=energies["ecore"],
        last_stage_energy_change=energy_change,
        last_stage_max_occupation_change=change,
        six_orbital_subspace_last_stage_max_sine_angle=subspace_change,
        error_estimate_caveat="Successive solver stages provide convergence evidence, not a rigorous error bound",
        fcidump_sha256=energies["fcidump_sha256"],
    )
    (path / "density_analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["provisional", "production"], required=True)
    args = ap.parse_args()
    ref = analyze(ROOT / "reference_zhendong")
    reconstructed = analyze(
        ROOT / ("cas20_provisional" if args.stage == "provisional" else "cas20")
    )
    difference = np.array(reconstructed["occupations"]) - ref["occupations"]
    reference_core = -4976.26532397
    comparison = dict(
        stage=args.stage,
        reference=ref,
        reconstructed=reconstructed,
        occupation_differences_reconstructed_minus_reference=difference.tolist(),
        maximum_absolute_difference=float(abs(difference).max()),
        rms_difference=float(np.sqrt(np.mean(difference**2))),
        reference_external_core=reference_core,
        reference_total_with_external_core=ref["energy_total"] + reference_core,
        total_energy_difference=reconstructed["energy_total"]
        - (ref["energy_total"] + reference_core),
        interpretation="Chemical and many-body diagnostics; neither energy nor occupation agreement establishes identical Hamiltonians or AO subspaces",
        selection_driven_by_energy_matching=False,
    )
    out = ROOT / "validation"
    (out / f"occupation_comparison_{args.stage}.json").write_text(
        json.dumps(comparison, indent=2) + "\n"
    )
    with (out / f"occupation_comparison_{args.stage}.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(["descending_index", "reference", "reconstructed", "difference"])
        for i, (a, b, d) in enumerate(
            zip(ref["occupations"], reconstructed["occupations"], difference)
        ):
            w.writerow([i, a, b, d])
    if args.stage == "production":
        subprocess.run(
            [sys.executable, str(CODE_ROOT / "11_assessment/sector_probability.py")],
            check=True,
        )
        sector = json.loads((ROOT / "cas20/sector_probability.json").read_text())
        smaller = analyze(ROOT / "cas14")
        provisional = analyze(ROOT / "cas20_provisional")
        dm = np.load(ROOT / "cas20/rdm1.npy")
        occ, u = np.linalg.eigh((dm + dm.T) / 2)
        u = u[:, ::-1]
        singular = np.linalg.svd(u[:6, :6], compute_uv=False)
        stability = float(np.sqrt(max(0.0, 1 - singular.min() ** 2)))
        dm2 = np.load(ROOT / "cas20/rdm2.npy")
        double = np.array([dm2[i, i, i, i] / 2 for i in range(6)])
        diag = np.diag(dm)[:6]
        compression = dict(
            delta_energy_freeze=smaller["energy_total"] - reconstructed["energy_total"],
            actual_K6_double_sector_probability=sector["probability"],
            actual_sector_probability_artifact="cas20/sector_probability.json",
            selected_fixed_basis_occupations=diag.tolist(),
            single_orbital_double_occupation_probabilities=double.tolist(),
            sector_weight_lower_bound_from_density=max(
                0.0, float(1 - np.sum(2 - diag))
            ),
            sector_weight_lower_bound_from_double_occupation=max(
                0.0, float(1 - np.sum(1 - double))
            ),
            sector_bound_interpretation="Projection onto fixed K6-double sector of CAS20; not fidelity with reoptimized CAS14",
            sorted_occupation_change_from_provisional=float(
                abs(
                    np.array(reconstructed["occupations"]) - provisional["occupations"]
                ).max()
            ),
            K6_production_vs_chosen_subspace_max_sine_angle=stability,
            variational_order_consistent_with_two_1e_5_Eh_solver_budgets=bool(
                smaller["energy_total"] - reconstructed["energy_total"] >= -2e-5
            ),
            status="pending_scientific_assessment",
            acceptance_energy_cutoff=None,
        )
        assert (
            sector["probability"]
            >= compression["sector_weight_lower_bound_from_double_occupation"] - 1e-7
        )
        (out / "compression_assessment.json").write_text(
            json.dumps(compression, indent=2) + "\n"
        )
    print(
        json.dumps(
            {
                "stage": args.stage,
                "occupation_max_difference": comparison["maximum_absolute_difference"],
                "total_energy_difference": comparison["total_energy_difference"],
            }
        )
    )


if __name__ == "__main__":
    main()
