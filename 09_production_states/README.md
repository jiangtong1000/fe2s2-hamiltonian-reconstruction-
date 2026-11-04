# 09_production_states: Production CAS20 and CAS14 states

Solve the production singlet CAS20 and candidate CAS14 using the final Hamiltonians.

## Method and scripts

The shared `run_dmrg.py` driver runs both production FCIDUMPs. The default schedule reaches bond dimension 2000. Inspect convergence, then use `--action refine-cas20 --bonds 3000 4000 5000` to continue; `refine-cas14` is available separately. Refinement preserves a copy of the prior state directory, so reserve substantial scratch space.

Do not choose a final state based only on reported sweep tolerance. Compare successive energies and sorted occupations, check 1-/2-RDM traces and energy reconstruction, reload the saved MPS, and evaluate the provisional trial in the production Hamiltonian. The target is 1e-5 Eh for new correlated energies, tightening when occupations or differences require it.

Scripts: [run_dmrg.py](../_shared/run_dmrg.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 09_production_states/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 09_production_states/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.
- `refine-cas20`: default refinement bonds [3000, 4000, 5000].
- `refine-cas14`: default refinement bonds [2500].

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 16 CPU(s), 160 GiB RAM, walltime 12:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `cas20/FCIDUMP`
- `cas14/FCIDUMP`

Key generated outputs (other diagnostic/state files are also written):

- `cas20/energies.json`
- `cas14/energies.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

Recorded CAS20: −5092.872362559896 Eh; last bond 4000→5000 change −5.31367e-6 Eh; maximum occupation change 5.19654e-5. Recorded CAS14: −5092.838193351019 Eh. Convergence histories, occupations, and compressed RDMs are retained. Successive-bond differences are empirical evidence, not rigorous error bars.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/cas14/convergence.json](results/cas14/convergence.json)
- [results/cas14/energies.json](results/cas14/energies.json)
- [results/cas14/occupations.txt](results/cas14/occupations.txt)
- [results/cas14_rdms.npz](results/cas14_rdms.npz)
- [results/cas20/convergence.json](results/cas20/convergence.json)
- [results/cas20/energies.json](results/cas20/energies.json)
- [results/cas20/occupations.txt](results/cas20/occupations.txt)
- [results/cas20_rdms.npz](results/cas20_rdms.npz)
- [results/validation/cas14_convergence_review.json](results/validation/cas14_convergence_review.json)
- [results/validation/cas20_stage5_review.json](results/validation/cas20_stage5_review.json)

## Reuse for a similar task

Expect schedules and resources to vary. Preserve the actual state/RDM provenance when resuming or recovering densities; a new random seed used for RDM recovery does not identify the earlier optimization seed.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
