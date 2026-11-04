"""Summarize this run without automatically accepting a compression candidate."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_shared"))
from workflow_paths import DATA_ROOT
import json


def main():
    def read(name):
        return json.loads((DATA_ROOT / name).read_text())

    comparison = read("validation/occupation_comparison_production.json")
    compression = read("validation/compression_assessment.json")
    validation = read("validation/validation_metrics.json")
    fci = read("validation/cas14_fci_comparison.json")
    cross = read("validation/cas20_cross_basis_energy.json")
    result = {
        "status": "computed_diagnostics_require_scientific_review",
        "numerical_validation_status": validation["status"],
        "cas14_fci_check": fci,
        "reference_total_energy_difference_Eh": comparison["total_energy_difference"],
        "sorted_occupation_max_difference": comparison["maximum_absolute_difference"],
        "compression_energy_cost_Eh": compression["delta_energy_freeze"],
        "K6_double_sector_probability": compression[
            "actual_K6_double_sector_probability"
        ],
        "cross_basis_check": cross,
        "compression_automatically_accepted": False,
        "interpretation": "Energy and natural occupations are diagnostics, not proof of historical Hamiltonian identity. Sector weight is not fidelity with reoptimized CAS14. Review convergence, chemistry, and intended application.",
    }
    (DATA_ROOT / "validation/scientific_assessment.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
