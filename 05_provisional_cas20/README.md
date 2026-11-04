# 05_provisional_cas20: Direct-AO CAS20 and basis-defining density

Compute a correlated density in the initial localized CAS20 basis before choosing a compression basis.

## Method and scripts

`prepare_provisional.py` transforms AO integrals directly into the initial 20-orbital space, contracting all 78 doubly occupied orbitals outside that CAS into its scalar and one-electron operator. This provisional CAS construction does not define the full Hamiltonian’s frozen core: the full model later contracts only F44.

The shared DMRG driver solves this provisional singlet and produces the 1-RDM used in step 06. Review energy, occupation, and prospective six-orbital-subspace convergence. `--action refine --bonds 1600 2000` resumes the saved state when needed. Do not rotate/finalize the basis while this density is still changing materially.

Scripts: [prepare_provisional.py](prepare_provisional.py), [run_dmrg.py](../_shared/run_dmrg.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 05_provisional_cas20/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 05_provisional_cas20/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.
- `refine`: default refinement bonds [1600, 2000].

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 16 CPU(s), 96 GiB RAM, walltime 12:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `orbitals/chemical_review.json`

Key generated outputs (other diagnostic/state files are also written):

- `cas20_provisional/energies.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

The archived direct-AO provisional FCIDUMP, energies, convergence, compressed RDMs, and basis-convergence review are retained. The generating density is explicitly labeled provisional and remains distinct from the final production density.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/cas20_provisional/FCIDUMP](results/cas20_provisional/FCIDUMP)
- [results/cas20_provisional/convergence.json](results/cas20_provisional/convergence.json)
- [results/cas20_provisional/energies.json](results/cas20_provisional/energies.json)
- [results/cas20_provisional/hamiltonian_provenance.json](results/cas20_provisional/hamiltonian_provenance.json)
- [results/cas20_provisional_rdms.npz](results/cas20_provisional_rdms.npz)
- [results/validation/occupation_comparison_provisional.json](results/validation/occupation_comparison_provisional.json)
- [results/validation/production_basis_review.json](results/validation/production_basis_review.json)

## Reuse for a similar task

Keep the independent direct-AO construction route. An active-space RDM alone is insufficient to regenerate full occupied–virtual couplings.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
