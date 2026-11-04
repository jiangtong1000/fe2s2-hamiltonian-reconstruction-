# 10_independent_validation: Direct integrals, embedding and CAS14 FCI

Check integrals, frozen-core contractions, embedding, and external-virtual couplings by independent routes.

## Method and scripts

`run_fci14.py` solves CAS14 with an independent iterative FCI solver. `validate_full.py` checks basis orthonormality/transformations, nested maps and electron identities, CAS FCIDUMP round trips, direct AO→CAS integrals, sequential freezing, factor-threshold sensitivity, and full embedded RDM energies. It now returns a failure exit status when checks fail or required states are missing.

The full-space checks go beyond CAS reconstruction: independently computed occupied/active–V h1 blocks, nine sampled ERI block families, and three determinants mixing occupied orbitals with external virtuals. The saved arrays are reloaded for these checks. The ERI sampling uses a fixed seed and is not an exhaustive check of all four-index elements.

`validate_fci14.py` compares FCI and DMRG energies/occupations and embeds the FCI RDMs through the full Hamiltonian. `--action smoke` runs small contraction, H2 factorization, and projector tests.

These are independent transformations/contractions and, for CAS14, solvers. They share the same PySCF AO integrals, geometry, basis, and physical hcore; they are not an independent electronic-structure-package calculation. They validate our model internally, not identity with the historical Hamiltonian.

Scripts: [ao_cholesky.py](../07_ao_cholesky/ao_cholesky.py), [run_fci14.py](run_fci14.py), [test_contractions.py](test_contractions.py), [validate_fci14.py](validate_fci14.py), [validate_full.py](validate_full.py), [sector_probability.py](../11_assessment/sector_probability.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 10_independent_validation/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 10_independent_validation/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.
- `smoke`.

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 16 CPU(s), 64 GiB RAM, walltime 04:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `cas20/rdm2.npy`
- `cas14/rdm2.npy`
- `reference_zhendong/rdm1.npy`

Key generated outputs (other diagnostic/state files are also written):

- `validation/validation_metrics.json`
- `validation/cas14_fci_comparison.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

Recorded checks passed. Representative errors: CAS20 h1 ~4.97e-10 Eh, ERI ~6.2e-11 Eh, scalar ~4.64e-9 Eh; sequential freezing h1 ~1.69e-14 Eh and scalar ~3.64e-12 Eh. CAS20 full embedding differs from its solver energy by 9.09e-12 Eh; direct-AO expectation differs by about 1.38e-8 Eh. Virtual-mixed determinant discrepancies are below 1.9e-8 Eh. CAS14 independent FCI is −5092.8381933518895 Eh, ~8.7e-10 Eh from DMRG.

Recorded final validation reused unchanged baseline integral checks only after verifying input hashes, while reevaluating the final states. This portable step recomputes the full validation rather than assuming those old checks apply.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/cas14/fci_validation/energies.json](results/cas14/fci_validation/energies.json)
- [results/validation/cas14_fci_comparison.json](results/validation/cas14_fci_comparison.json)
- [results/validation/report.txt](results/validation/report.txt)
- [results/validation/validation_metrics.json](results/validation/validation_metrics.json)

## Reuse for a similar task

Keep tests for the external space as well as active projections. Investigate failed checks rather than loosening thresholds to accommodate them.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
