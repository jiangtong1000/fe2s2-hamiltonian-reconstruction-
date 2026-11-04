# 04_reference_dmrg: Historical Hamiltonian: singlet DMRG and RDMs

Establish the historical singlet reference using the same solver and RDM conventions as our reconstruction.

## Method and scripts

`run.py` calls `../_shared/run_dmrg.py` with the historical FCIDUMP directory as its working directory. It uses block2 SU2 singlet symmetry, seed 10405, input orbital order, a two-site optimization followed by one-site finishing sweeps, and saves native GS states, convergence history, and spin-summed 1-/2-RDMs.

The corrected shared driver restores the two-site representation at every bond stage. A fresh run is a reproducible construction recipe, not a byte-for-byte replay of earlier debugging and restart history. Historical submitted snapshots are retained in `../provenance/submitted/`.

Use `--action refine --bonds 1600 2000` only after examining density convergence; the initial schedule already reaches 2000.

Scripts: [run_dmrg.py](../_shared/run_dmrg.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 04_reference_dmrg/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 04_reference_dmrg/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.
- `refine`: default refinement bonds [1600, 2000].

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 16 CPU(s), 96 GiB RAM, walltime 12:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `reference_zhendong/FCIDUMP`

Key generated outputs (other diagnostic/state files are also written):

- `reference_zhendong/energies.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

Recorded matched-reference active energy: −116.605608594803 Eh. Adding the external core once gives −5092.870932564803 Eh. The maximum occupation change at the last recorded refinement is approximately 1.55e-5. Compressed RDMs and convergence history are committed; native MPS states remain external.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/reference_zhendong/convergence.json](results/reference_zhendong/convergence.json)
- [results/reference_zhendong/energies.json](results/reference_zhendong/energies.json)
- [results/reference_zhendong/occupations.txt](results/reference_zhendong/occupations.txt)
- [results/reference_zhendong/provenance.json](results/reference_zhendong/provenance.json)
- [results/reference_zhendong_rdms.npz](results/reference_zhendong_rdms.npz)
- [results/validation/reference_convergence_review.json](results/validation/reference_convergence_review.json)

## Reuse for a similar task

Preserve the reference FCIDUMP scalar convention, input ordering, state symmetry, and RDM convention. Matching solver settings helps interpret differences but does not establish equivalent active spaces.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
