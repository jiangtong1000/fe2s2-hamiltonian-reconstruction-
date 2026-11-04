"""Pivoted modified Cholesky of AO Coulomb integrals without full AO ERIs.

Pairs are unweighted lower-triangular AO pairs; factors reconstruct (mu nu|la si)
by an ordinary dot product. No svec sqrt(2) convention is used.
"""

from pathlib import Path as _Path
import sys as _sys

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "_shared"))
from workflow_paths import DATA_ROOT
import argparse
import json
import os
import time
import numpy as np
from pyscf import gto, lib


def decompose(mol, threshold=1e-9, max_rank=None, progress=None):
    n = mol.nao_nr()
    rows, cols = np.tril_indices(n)
    npair = len(rows)
    max_rank = max_rank or 15 * n
    # Allocation is O(nAO^2 * rank), not a four-index molecular ERI tensor.
    factors = np.empty((max_rank, npair), dtype=np.float64)
    diag = np.empty((n, n))
    offsets = mol.ao_loc_nr()
    for sh in range(mol.nbas):
        a, b = offsets[sh : sh + 2]
        buf = mol.intor(
            "int2e_sph",
            shls_slice=(sh, sh + 1, 0, mol.nbas, sh, sh + 1, 0, mol.nbas),
            aosym="s1",
        )
        diag[a:b] = buf.reshape((b - a) * n, (b - a) * n).diagonal().reshape(b - a, n)
    residual = diag[rows, cols].copy()
    milestones = {}
    start = time.time()
    rank = 0
    while True:
        pivot = int(residual.argmax())
        error = float(residual[pivot])
        for tol in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]:
            if error <= tol and str(tol) not in milestones:
                milestones[str(tol)] = rank
        if error <= threshold:
            break
        if rank >= max_rank:
            raise RuntimeError("Cholesky rank capacity exceeded; increase max_rank")
        mu, nu = int(rows[pivot]), int(cols[pivot])
        smu = int(np.searchsorted(offsets, mu, side="right") - 1)
        snu = int(np.searchsorted(offsets, nu, side="right") - 1)
        buf = mol.intor(
            "int2e_sph",
            shls_slice=(0, mol.nbas, 0, mol.nbas, smu, smu + 1, snu, snu + 1),
            aosym="s1",
        )
        column = buf[:, :, mu - offsets[smu], nu - offsets[snu]][rows, cols].copy()
        if rank:
            column -= factors[:rank, pivot] @ factors[:rank]
        factors[rank] = column / np.sqrt(error)
        residual -= factors[rank] ** 2
        rank += 1
        if rank % 25 == 0 and progress:
            progress(
                dict(
                    rank=rank,
                    pivot_residual=error,
                    residual_max=float(residual.max()),
                    residual_min=float(residual.min()),
                    elapsed_seconds=time.time() - start,
                )
            )
        if residual.min() < -1e-7:
            raise RuntimeError(f"Unexpected negative Coulomb residual {residual.min()}")
    return factors[:rank], dict(
        rank=rank,
        nao=n,
        npair=npair,
        threshold=threshold,
        residual_max=float(residual.max()),
        residual_min=float(residual.min()),
        ranks_at_threshold=milestones,
        elapsed_seconds=time.time() - start,
        pair_convention="numpy.tril_indices; unweighted symmetric AO pairs",
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--threshold", type=float, default=1e-9)
    p.add_argument("--smoke", action="store_true")
    args = p.parse_args()
    lib.num_threads(int(os.environ.get("OMP_NUM_THREADS", "1")))
    root = DATA_ROOT
    if args.smoke:
        mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="cc-pvdz", verbose=0)
        l, meta = decompose(mol, 1e-10, max_rank=mol.nao_nr() ** 2)
        eri = mol.intor("int2e", aosym="s4")
        err = float(abs(eri - l.T @ l).max())
        assert err < 1.01e-10, err
        meta["full_eri_reconstruction_max_error"] = err
        print(json.dumps(meta, indent=2))
        return
    out = root / "ao_factorization"
    out.mkdir(exist_ok=True)
    mol = gto.M(
        atom=str(root / "geometry.xyz"),
        basis=json.loads((root / "inputs/basis.json").read_text()),
        charge=-2,
        spin=10,
        unit="Angstrom",
        cart=False,
        verbose=4,
        max_memory=48000,
    )

    def progress(meta):
        print(json.dumps(meta), flush=True)
        (out / "progress.json").write_text(json.dumps(meta, indent=2) + "\n")

    factors, meta = decompose(mol, args.threshold, progress=progress)
    np.save(out / "chol_ao_packed.npy", factors)
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print("AO_CHOLESKY_COMPLETE", json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
