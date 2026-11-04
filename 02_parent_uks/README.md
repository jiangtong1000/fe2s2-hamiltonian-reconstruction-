# 02_parent_uks: High-spin BP86 and parent natural orbitals

Build the high-spin molecular reference and a complete orthonormal parent natural-orbital basis.

## Method and scripts

`run_parent.py` runs unrestricted BP86/TZP-DKH with spin-free X2C, grid level 4, SCF energy tolerance 1e-10 Eh, and gradient tolerance 1e-6. The parent has S=5, multiplicity 11, PySCF spin=10, and (Nα,Nβ)=(98,88). This is scalar relativity, not spin–orbit coupling. The physical hcore does not include the DFT exchange-correlation potential.

Diagonalize S^(1/2)(Dα+Dβ)S^(1/2), then transform eigenvectors back with S^(-1/2). Save both AO spin densities, all coefficients, occupations, overlap, hcore, occupied UKS orbitals, and AO labels. These natural orbitals come from the **spin-summed high-spin UKS 1-RDM**, not a correlated singlet density. `df_total_energy` in the historical metadata is the DFT energy, not a density-fitting flag.

The `newton` action chooses an alternative fresh second-order SCF in an unused parent output directory. It is not an automatic fallback or a claim about the archived SCF algorithm.

Scripts: [run_parent.py](run_parent.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 02_parent_uks/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 02_parent_uks/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.
- `newton`.

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 8 CPU(s), 64 GiB RAM, walltime 02:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `geometry.xyz`
- `inputs/basis.json`

Key generated outputs (other diagnostic/state files are also written):

- `parent/scf.chk`
- `parent/ao_metadata.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

662 spherical AOs; 186 electrons; complete parent basis retained. The recorded parent metadata and occupation table are in `results/parent/`. The checkpoint and parent arrays are omitted from Git and rebuilt by this step.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/parent/ao_metadata.json](results/parent/ao_metadata.json)
- [results/parent/parent_natural_occupations.txt](results/parent/parent_natural_occupations.txt)

## Reuse for a similar task

Review relativity, AO linear dependence, spin state, SCF solution, and physical hcore treatment for the new molecule. Never compare a DFT total energy with a wavefunction Hamiltonian expectation as if they were identical.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
