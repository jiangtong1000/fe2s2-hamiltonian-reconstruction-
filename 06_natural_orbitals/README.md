# 06_natural_orbitals: Correlated natural orbitals and nested partition

Define the correlated natural-orbital rotation and preserve exact nesting and trial semantics.

## Method and scripts

`finalize_orbitals.py` diagonalizes the provisional **singlet CAS(30e,20o) DMRG 1-RDM**. It sorts occupations descending, fixes phases by the largest AO coefficient, and resolves numerical occupation degeneracies using projected physical hcore. These orbitals are different from the parent UKS natural orbitals.

The six most occupied orbitals form a candidate K6; the remaining fourteen form A14. The initial CAS20 subspace is only rotated. The complete retained basis is [I34 | K6 | A14 | V564], with F44 separately recorded. `describe_final_basis.py` writes shell characters and both trial embeddings, including index, spin-orbital, and determinant-phase conventions.

K6 and I are fixed occupied in the appropriate trials; they remain in the full propagation space. V is empty in those trials, but retained with its full interactions. No CASSCF/DMRG-SCF orbital optimization is introduced.

Scripts: [describe_final_basis.py](describe_final_basis.py), [finalize_orbitals.py](finalize_orbitals.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 06_natural_orbitals/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 06_natural_orbitals/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 4 CPU(s), 8 GiB RAM, walltime 01:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `cas20_provisional/rdm1.npy`

Key generated outputs (other diagnostic/state files are also written):

- `orbitals/index_maps.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

Exact orbital nesting A14 ⊂ A20 ⊂ H_full; 618 retained orbitals and 98 electrons. `results/final_basis.npz` retains final full/frozen/initial-CAS coefficients and the CAS20 rotation. Index maps, K6 chemical characters, and embedding semantics are also committed. The fixed rotation was generated from the provisional density; later production occupations validate its stability.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/final_basis.npz](results/final_basis.npz)
- [results/orbitals/K6_character.json](results/orbitals/K6_character.json)
- [results/orbitals/index_maps.json](results/orbitals/index_maps.json)
- [results/orbitals/natural_orbital_diagnostics.json](results/orbitals/natural_orbital_diagnostics.json)
- [results/orbitals/trial_embeddings.json](results/orbitals/trial_embeddings.json)

## Reuse for a similar task

Inspect the occupation boundary, subspace stability, chemical character, and eventual energy cost. A high occupation threshold is a screen, not automatic permission to freeze six orbitals.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
