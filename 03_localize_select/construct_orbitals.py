"""Chemical selection in a complete, orthonormal split-localized UNO basis.

No energy-based selection. Deep core is projected inside the near-double UNO
block; occupied valence and singly occupied blocks are localized separately.
"""

from pathlib import Path as _Path
import sys as _sys

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "_shared"))
from workflow_paths import DATA_ROOT
import argparse
import csv
import json
import os
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment
from pyscf import scf, lo, lib

ROOT = DATA_ROOT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nfrozen", type=int, choices=[0, 44], required=True)
    args = ap.parse_args()
    lib.num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))
    out = ROOT / "orbitals"
    out.mkdir(exist_ok=True)
    mol = scf.chkfile.load_mol(str(ROOT / "parent/scf.chk"))
    c = np.load(ROOT / "parent/mo_coeff.npy")
    no = np.load(ROOT / "parent/parent_natural_occupations.npy")
    s = np.load(ROOT / "parent/overlap_ao.npy")
    h = np.load(ROOT / "parent/hcore_ao.npy")
    labels = mol.ao_labels(fmt=False)
    meta = lo.orth.orth_ao(mol, "meta_lowdin", s=s)
    groups = {
        "double": np.flatnonzero(no > 1.5),
        "single": np.flatnonzero((no >= 0.5) & (no <= 1.5)),
        "virtual": np.flatnonzero(no < 0.5),
    }
    counts = {key: len(val) for key, val in groups.items()}
    assert counts == {"double": 88, "single": 10, "virtual": 564}, counts
    # Explicit deep-core model: Fe/S 1s,2s,2p and C 1s. Fe3s/3p remain valence.
    core_aos = [
        i
        for i, (atom, el, shell, component) in enumerate(labels)
        if (el in ["Fe", "S"] and shell in ["1s", "2s", "2p"])
        or (el == "C" and shell == "1s")
    ]
    assert len(core_aos) == 44, len(core_aos)
    cd = c[:, groups["double"]]
    proj = meta[:, core_aos].T @ s @ cd
    core_eig, u = eigh(proj.T @ proj)
    u = u[:, ::-1]
    core_eig = core_eig[::-1]
    core0 = cd @ u[:, :44]
    val0 = cd @ u[:, 44:]

    def pm(block, label):
        loc = lo.PM(mol, block, pop_method="mulliken")
        loc.exponent = 2
        loc.conv_tol = 1e-10
        loc.conv_tol_grad = 1e-6
        loc.max_cycle = 150
        loc.init_guess = "atomic"
        answer = loc.kernel()
        np.save(out / (label + "_localized.npy"), answer)
        return answer, dict(
            pop_method="mulliken",
            exponent=2,
            conv_tol=1e-10,
            conv_tol_grad=1e-6,
            max_cycle=150,
            initial_guess="atomic",
            final_objective=float(loc.cost_function()),
            final_gradient_norm=float(np.linalg.norm(loc.get_grad())),
        )

    core, settings_core = pm(core0, "deep_core")
    val, settings_val = pm(val0, "occupied_valence")
    single, settings_single = pm(c[:, groups["single"]], "single")
    # Meta-Lowdin weights provide separate reproducible character diagnostics.
    pv = meta.T @ s @ val
    ps = meta.T @ s @ single

    def ao_idx(atom, shell):
        return [i for i, x in enumerate(labels) if x[0] == atom and x[2] == shell]

    geo = json.loads((ROOT / "sources/geometry_provenance.json").read_text())
    fes = geo["fe_indices_zero_based"]
    bridges = geo["bridge_s_indices_zero_based"]
    terminals = geo["terminal_s_indices_zero_based"]
    xyz = mol.atom_coords(unit="Angstrom")

    def sweight(atom, directional=False):
        ids = ao_idx(atom, "3p")
        assert len(ids) == 3, (atom, ids)
        if directional:
            nearest = min(fes, key=lambda fe: np.linalg.norm(xyz[fe] - xyz[atom]))
            direction = xyz[nearest] - xyz[atom]
            direction /= np.linalg.norm(direction)
            unit = {"x": direction[0], "y": direction[1], "z": direction[2]}
            coeff = np.array([unit[labels[k][3]] for k in ids])
            return (coeff @ pv[ids]) ** 2
        return np.sum(pv[ids] ** 2, axis=0)

    # Six bridging orbitals (three per S), four terminal Fe-directed orbitals.
    slots = [(atom, False) for atom in bridges for _ in range(3)] + [
        (atom, True) for atom in terminals
    ]
    scores = np.array([sweight(atom, directed) for atom, directed in slots])
    slot_ids, chosen = linear_sum_assignment(-scores)
    assert np.array_equal(slot_ids, np.arange(10))
    fe_scores = np.array(
        [np.sum(ps[ao_idx(atom, "3d")] ** 2, axis=0) for atom in fes for _ in range(5)]
    )
    _, chosen_d = linear_sum_assignment(-fe_scores)
    s_for_atom = {
        atom: list(chosen[:6][np.repeat(np.array(bridges), 3) == atom])
        for atom in bridges
    }
    term_choice = dict(zip(terminals, chosen[6:]))
    left = [
        atom
        for atom in terminals
        if np.linalg.norm(xyz[atom] - xyz[fes[0]])
        < np.linalg.norm(xyz[atom] - xyz[fes[1]])
    ]
    right = [atom for atom in terminals if atom not in left]
    assert len(left) == len(right) == 2
    a_cols = [val[:, term_choice[t]] for t in left]
    a_cols += [single[:, j] for j in chosen_d[:5]]
    a_cols += [val[:, s_for_atom[b][j]] for j in range(3) for b in bridges]
    a_cols += [single[:, j] for j in chosen_d[5:]]
    a_cols += [val[:, term_choice[t]] for t in right]
    ca = np.column_stack(a_cols)
    inactive_val = val[:, [j for j in range(val.shape[1]) if j not in chosen]]
    frozen = core if args.nfrozen else core[:, :0]
    inactive = inactive_val if args.nfrozen else np.column_stack([core, inactive_val])
    virtual = c[:, groups["virtual"]]
    complete = np.column_stack([frozen, inactive, ca, virtual])
    err = float(abs(complete.T @ s @ complete - np.eye(len(no))).max())
    assert err < 1e-10, err
    assert args.nfrozen + inactive.shape[1] == 78
    for name, a in dict(
        mo_coeff_partition_initial=complete,
        mo_coeff_frozen_core=frozen,
        mo_coeff_20_initial=ca,
        selection_from_parent=c.T @ s @ complete,
        character_meta_lowdin=meta.T @ s @ complete,
    ).items():
        np.save(out / (name + ".npy"), a)
    np.savez(
        out / "selection_diagnostics.npz",
        core_projector_eigenvalues=core_eig,
        ligand_selection_scores=scores,
        ligand_selected_columns=chosen,
        fe_selection_scores=fe_scores,
        fe_selected_columns=chosen_d,
        parent_natural_occupations=no,
    )
    meta_pops = (meta.T @ s @ complete) ** 2
    atom_pops = np.array(
        [
            meta_pops[[i for i, x in enumerate(labels) if x[0] == a]].sum(axis=0)
            for a in range(mol.natm)
        ]
    )
    dparent = np.load(ROOT / "parent/density_alpha_ao.npy") + np.load(
        ROOT / "parent/density_beta_ao.npy"
    )
    pops = np.diag(complete.T @ s @ dparent @ s @ complete)
    e_h = np.diag(complete.T @ h @ complete)
    fields = [
        "partition_index_zero_based",
        "block",
        "parent_UKS_occupation",
        "physical_hcore_expectation_Eh",
        "center_atom_zero_based",
    ] + [f"atom{a}_{mol.atom_symbol(a)}_meta_lowdin_weight" for a in range(mol.natm)]
    with (out / "orbital_table_initial.csv").open("w") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for j in range(complete.shape[1]):
            block = (
                "F" if j < args.nfrozen else "I" if j < 78 else "A20" if j < 98 else "V"
            )
            w.writerow(
                [
                    j,
                    block,
                    pops[j],
                    e_h[j],
                    int(atom_pops[:, j].argmax()),
                    *atom_pops[:, j],
                ]
            )
    report = dict(
        status="constructed_pending_chemical_review",
        nfrozen=args.nfrozen,
        ninactive=inactive.shape[1],
        nactive=20,
        nvirtual=virtual.shape[1],
        orthonormality_max_error=err,
        uno_blocks={k: v.tolist() for k, v in groups.items()},
        uno_block_thresholds=[0.5, 1.5],
        core_definition="Fe/S [Ne] and C1s"
        if args.nfrozen
        else "none; all 186 electrons retained",
        core_projector_smallest_selected=float(core_eig[43]),
        core_projector_largest_excluded=float(core_eig[44]),
        pm_settings=dict(
            core=settings_core, occupied_valence=settings_val, single=settings_single
        ),
        character_method="meta-Lowdin orthogonal AOs with ANO preorthogonalization",
        ligand_assignments=[
            dict(
                atom=int(a),
                terminal_fe_directed=bool(t),
                valence_column=int(j),
                score=float(scores[k, j]),
            )
            for k, ((a, t), j) in enumerate(zip(slots, chosen))
        ],
        fe_assignments=[
            dict(
                atom=int(fes[k // 5]),
                single_column=int(j),
                score=float(fe_scores[k, j]),
            )
            for k, j in enumerate(chosen_d)
        ],
        production_basis_finalized=False,
    )
    (out / "selection_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print("SELECTION_READY_FOR_REVIEW", json.dumps(report))


if __name__ == "__main__":
    main()
