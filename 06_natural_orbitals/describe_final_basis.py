"""Character and embedding metadata for the final retained basis."""

from pathlib import Path as _Path
import sys as _sys

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "_shared"))
from workflow_paths import DATA_ROOT
import csv
import json
import numpy as np
from pyscf import scf

ROOT = DATA_ROOT


def main():
    out = ROOT / "orbitals"
    mol = scf.chkfile.load_mol(str(ROOT / "parent/scf.chk"))
    labels = mol.ao_labels(fmt=False)
    p = np.load(out / "character_meta_lowdin.npy")[:, 44:].copy()
    p[:, 34:54] = p[:, 34:54] @ np.load(out / "U_natural_20.npy")
    maps = json.loads((out / "index_maps.json").read_text())
    groups = [
        (a, shell)
        for a in range(mol.natm)
        for shell in dict.fromkeys(x[2] for x in labels if x[0] == a)
    ]
    weights = {
        f"atom{a}_{mol.atom_symbol(a)}_{shell}": np.sum(
            p[[i for i, x in enumerate(labels) if x[0] == a and x[2] == shell]] ** 2,
            axis=0,
        )
        for a, shell in groups
    }
    occ = json.loads((out / "natural_orbital_diagnostics.json").read_text())[
        "occupations"
    ]
    with (out / "orbital_character_final.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "retained_index_zero_based",
                "block",
                "provisional_CAS20_natural_occupation",
                *weights,
            ]
        )
        for i in range(p.shape[1]):
            block = next(b for b in ["I", "K6", "A14", "V"] if i in maps[b])
            w.writerow(
                [
                    i,
                    block,
                    occ[i - 34] if 34 <= i < 54 else "",
                    *[x[i] for x in weights.values()],
                ]
            )
    summary = []
    for i in maps["K6"]:
        characters = sorted(
            [(key, float(val[i])) for key, val in weights.items()],
            key=lambda item: -item[1],
        )
        summary.append(
            dict(
                retained_index=i,
                CAS20_index=i - 34,
                occupation=occ[i - 34],
                leading_shell_weights=characters[:8],
            )
        )
    (out / "K6_character.json").write_text(json.dumps(summary, indent=2) + "\n")
    embedding = {}
    for name, active, fixed in [
        ("cas20", maps["cas20_local_to_full_zero_based"], maps["I"]),
        ("cas14", maps["A14"], maps["I"] + maps["K6"]),
    ]:
        embedding[name] = dict(
            active_local_to_retained=active,
            fixed_double_retained=fixed,
            empty_retained=maps["V"],
            nalpha=49,
            nbeta=49,
            native_state_directory=f"{name}/native_state",
            native_state_tag="GS",
            state_symmetry="block2 SU2 singlet",
            native_state_requires_block2=True,
            spin_orbital_order=maps["spin_orbital_order"],
            determinant_embedding_phase_rule="Map active determinant creation operators, prepend fixed alpha then fixed beta creation operators, and multiply by parity of permutation to globally ascending retained spin-orbital order. A common overall phase is immaterial.",
            expectation_route="10_independent_validation/validate_full.py via _shared/integrals.py (repository paths): full-density Wick terms plus active 2-RDM cumulant",
            occupations_are_trial_properties_only=True,
            frozen_F_already_contracted=True,
        )
    (out / "trial_embeddings.json").write_text(json.dumps(embedding, indent=2) + "\n")
    print("FINAL_BASIS_CHARACTER_AND_EMBEDDINGS_WRITTEN")


if __name__ == "__main__":
    main()
