# 03_localize_select: Frozen core, split localization and chemical CAS20

Select a chemically defined CAS(30e,20o), retaining a complete complement and identifying genuine frozen core.

## Method and scripts

`construct_orbitals.py --nfrozen 44` partitions parent occupations into double (>1.5), single (0.5–1.5), and virtual (<0.5) blocks. It identifies the deep-core projector inside the double block and split-localizes valence blocks with Pipek–Mezey. The implementation uses Mulliken PM populations, exponent 2, energy/objective tolerance 1e-10, and gradient tolerance 1e-6. These are reproducible reconstruction choices; exact historical settings are unverified.

`review_selection.py` checks meta-Löwdin shell weights and bonded-ligand directionality. The selected space is ten Fe 3d-derived orbitals, six bridging S 3p-derived orbitals, and four terminal S 3p orbitals directed toward Fe. Molecular ligand mixing is expected. The initial ordering follows terminal2 / Fe1(5) / bridge6 / Fe2(5) / terminal2. No orbital is selected by agreement with the historical DMRG energy.

Freeze only Fe/S [Ne] and C 1s: 44 spatial orbitals. Fe 3s/3p remain correlated in the full Hamiltonian. Inspect the chemical review before continuing; do not bypass failed assertions by adjusting energy indices.

Scripts: [construct_orbitals.py](construct_orbitals.py), [review_selection.py](review_selection.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 03_localize_select/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 03_localize_select/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 8 CPU(s), 32 GiB RAM, walltime 02:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `parent/scf.chk`

Key generated outputs (other diagnostic/state files are also written):

- `orbitals/chemical_review.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

Chemical review passed. The initial complete partition is F44 / I34 / A20 / V564. `results/CAS20_shell_characters.csv` retains the active rows of the full shell-character table; selection and review JSON preserve the criteria and conclusions.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/CAS20_shell_characters.csv](results/CAS20_shell_characters.csv)
- [results/orbitals/chemical_review.json](results/orbitals/chemical_review.json)
- [results/orbitals/initial_cas20_roles.json](results/orbitals/initial_cas20_roles.json)
- [results/orbitals/selection_report.json](results/orbitals/selection_report.json)
- [results/validation/initial_transformation_roundtrip.json](results/validation/initial_transformation_roundtrip.json)

## Reuse for a similar task

Adapt atom indices, shell projectors, terminal bond vectors, block counts, and assignment rules explicitly. Chemical active-space selection must be reviewed anew for a different molecule.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
