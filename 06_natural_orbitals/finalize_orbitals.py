"""CAS20 correlated natural orbitals; full complement preserved exactly."""

from pathlib import Path as _Path
import sys as _sys

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "_shared"))
from workflow_paths import DATA_ROOT
import hashlib
import json
import os
import numpy as np
from scipy.linalg import eigh

ROOT = DATA_ROOT


def main():
    out = ROOT / "orbitals"
    selection = json.loads((out / "selection_report.json").read_text())
    assert selection["status"] == "chemical_review_passed"
    c = np.load(out / "mo_coeff_partition_initial.npy")
    s = np.load(ROOT / "parent/overlap_ao.npy")
    rdm_path = ROOT / "cas20_provisional/rdm1.npy"
    dm = np.load(rdm_path)
    assert abs(dm - dm.T).max() < 1e-8
    assert abs(np.trace(dm) - 30) < 1e-7
    occ, u = eigh((dm + dm.T) / 2)
    occ = occ[::-1]
    u = u[:, ::-1]
    assert occ.min() > -1e-8 and occ.max() < 2 + 1e-8
    # Within numerical degeneracies, diagonalize the saved physical hcore projection.
    ca = c[:, 78:98]
    h = ca.T @ np.load(ROOT / "parent/hcore_ao.npy") @ ca
    start = 0
    degeneracies = []
    while start < 20:
        end = start + 1
        while end < 20 and abs(occ[end] - occ[start]) < 1e-10:
            end += 1
        if end - start > 1:
            _, v = eigh(u[:, start:end].T @ h @ u[:, start:end])
            u[:, start:end] = u[:, start:end] @ v
            degeneracies.append(list(range(start, end)))
        start = end
    ca_no = ca @ u
    phase = np.where(ca_no[np.argmax(abs(ca_no), axis=0), np.arange(20)] < 0, -1.0, 1.0)
    u *= phase
    ca_no *= phase
    nf = selection["nfrozen"]
    ni = selection["ninactive"]
    f = c[:, :nf]
    inactive = c[:, nf:78]
    v = c[:, 98:]
    full = np.column_stack([inactive, ca_no, v])
    err = float(abs(full.T @ s @ full - np.eye(full.shape[1])).max())
    assert err < 1e-10, err
    for name, a in dict(
        U_natural_20=u,
        mo_coeff_20_final=ca_no,
        mo_coeff_6=ca_no[:, :6],
        mo_coeff_14=ca_no[:, 6:],
        mo_coeff_full_final=full,
        dmrg_rdm1_20=dm,
        mo_coeff_frozen_core=f,
    ).items():
        np.save(out / (name + ".npy"), a)
    np.savez(
        out / "transformations.npz",
        selection_from_parent=np.load(out / "selection_from_parent.npy"),
        natural_rotation20=u,
        final_partition_from_initial=np.block(
            [
                [np.eye(78), np.zeros((78, 20)), np.zeros((78, c.shape[1] - 98))],
                [np.zeros((20, 78)), u, np.zeros((20, c.shape[1] - 98))],
                [np.zeros((c.shape[1] - 98, 98)), np.eye(c.shape[1] - 98)],
            ]
        ),
    )
    maps = dict(
        coordinate_system="retained_spatial_array_zero_based",
        F_index_space="initial complete localized partition [F|I|A20|V], not original parent UNO indices",
        cas20_rotation_convention="C20_final = C20_initial @ U_natural_20; this is a dense rotation, not a permutation",
        complete_basis_transform_convention="C_final_complete = C_parent_UNO @ selection_from_parent @ final_partition_from_initial (both transformations in transformations.npz)",
        complete_final_order=["F", "I", "K6", "A14", "V"],
        retained_complete_partition_slice_zero_based=[nf, c.shape[1]],
        F_parent_partition_zero_based=list(range(nf)),
        I=list(range(ni)),
        K6=list(range(ni, ni + 6)),
        A14=list(range(ni + 6, ni + 20)),
        V=list(range(ni + 20, full.shape[1])),
        cas20_local_to_full_zero_based=list(range(ni, ni + 20)),
        cas14_local_to_full_zero_based=list(range(ni + 6, ni + 20)),
        cas20_local_to_full_fcidump_one_based=list(range(ni + 1, ni + 21)),
        cas14_local_to_full_fcidump_one_based=list(range(ni + 7, ni + 21)),
        nelec_full=186 - 2 * nf,
        nelec_cas20=30,
        nelec_cas14=18,
        spin_orbital_order="all retained alpha ascending, then all retained beta ascending",
        determinant_phase_convention="creation operators ascending in declared spin-orbital order",
        fixed_trial_occupations_restrict_full_propagation=False,
    )
    (out / "index_maps.json").write_text(json.dumps(maps, indent=2) + "\n")
    diagnostic = dict(
        source="singlet CAS(30e,20o) DMRG 1-RDM",
        source_stage="provisional_direct_AO",
        generating_slurm_job_id=os.environ.get("SLURM_JOB_ID"),
        source_rdm_sha256=hashlib.sha256(rdm_path.read_bytes()).hexdigest(),
        occupation_sort="descending",
        occupations=occ.tolist(),
        selected_six=occ[:6].tolist(),
        delta6=float(np.sum(2 - occ[:6])),
        sector_weight_lower_bound=max(0.0, float(1 - np.sum(2 - occ[:6]))),
        gap_at_boundary=float(occ[5] - occ[6]),
        all_six_above_1p97=bool(np.all(occ[:6] >= 1.97)),
        degenerate_blocks=degeneracies,
        degeneracy_tolerance=1e-10,
        degeneracy_rule="diagonalize physical hcore within nearly degenerate occupation blocks",
        phase_rule="largest magnitude AO coefficient positive",
        orthonormality_max_error=err,
        offdiagonal_rdm_max=float(
            abs(u.T @ dm @ u - np.diag(np.diag(u.T @ dm @ u))).max()
        ),
        compression_status="candidate_pending_chemical_and_energy_assessment",
    )
    (out / "natural_orbital_diagnostics.json").write_text(
        json.dumps(diagnostic, indent=2) + "\n"
    )
    print("FINAL_BASIS_CANDIDATE", json.dumps(diagnostic))


if __name__ == "__main__":
    main()
