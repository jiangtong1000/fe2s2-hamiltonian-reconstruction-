# 11_assessment: Historical comparison and compression assessment

Assess chemical/many-body similarity and quantify the six-orbital compression without selecting by energy.

## Method and scripts

`analyze_states.py` applies the same sorted-occupation analysis to historical and reconstructed singlet RDMs, verifies stage consistency, and computes energy differences using the historical external core exactly once. `sector_probability.py` evaluates the actual probability of the fixed K6-double sector from the saved CAS20 state. `_shared/cross_basis_energy.py` checks a saved provisional trial in the production Hamiltonian by rotating integrals and, independently, densities.

`plot_diagnostics.py` generates shareable PDF/PNG character and occupation plots. `summarize.py` assembles diagnostics for a new run, explicitly leaving application-specific compression acceptance for scientific review. It does not stamp a fresh run with the archived conclusion.

Scripts: [analyze_states.py](analyze_states.py), [plot_diagnostics.py](plot_diagnostics.py), [summarize.py](summarize.py), [cross_basis_energy.py](../_shared/cross_basis_energy.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 11_assessment/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 11_assessment/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 4 CPU(s), 16 GiB RAM, walltime 02:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `cas20/rdm2.npy`
- `cas14/rdm2.npy`
- `reference_zhendong/rdm1.npy`

Key generated outputs (other diagnostic/state files are also written):

- `validation/occupation_comparison_production.json`
- `validation/compression_assessment.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

Our recorded total energy is 1.429995 mEh below the matched historical reference. Maximum sorted occupation difference: 0.0019006599; RMS: 0.00090487646. Similar occupations and chemistry do not prove equal Hamiltonians or AO subspaces.

Fixing K6 costs 34.1692089 mEh. Its exact double-sector probability is 0.9403626731; this is not fidelity with the reoptimized CAS14 state. The six orbitals pass the 1.97 occupation screen, but CAS14 is an aggressive, quantified compression candidate—not an energy-equivalent replacement for CAS20. Application-specific VCE/AFQMC usefulness remains open.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/cas20/sector_probability.json](results/cas20/sector_probability.json)
- [results/figures/initial_CAS20_character.pdf](results/figures/initial_CAS20_character.pdf)
- [results/figures/occupation_comparison_production.pdf](results/figures/occupation_comparison_production.pdf)
- [results/figures/occupation_comparison_provisional.pdf](results/figures/occupation_comparison_provisional.pdf)
- [results/validation/cas20_cross_basis_energy.json](results/validation/cas20_cross_basis_energy.json)
- [results/validation/compression_assessment.json](results/validation/compression_assessment.json)
- [results/validation/occupation_comparison_production.csv](results/validation/occupation_comparison_production.csv)
- [results/validation/occupation_comparison_production.json](results/validation/occupation_comparison_production.json)
- [results/validation/scientific_assessment.json](results/validation/scientific_assessment.json)

## Reuse for a similar task

Compare chemically defined spaces without fitting orbitals to an energy target. Retain failed or costly compression outcomes as useful results.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
