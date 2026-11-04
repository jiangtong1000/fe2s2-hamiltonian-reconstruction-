# 07_ao_cholesky: Direct AO Cholesky factorization

Factor the AO Coulomb operator without storing the dense four-index tensor.

## Method and scripts

`ao_cholesky.py` performs pivoted modified Cholesky decomposition of spherical AO Coulomb integrals. The stopping residual is 1e-9 Eh. Lower-triangular AO pairs are unweighted (`numpy.tril_indices`): no sqrt(2) pair scaling is used. It records the rank at a looser 1e-8 threshold for later sensitivity checks.

`check_ao_factors.py` checks sampled reconstructed integrals against direct AO integrals and records the serialized array hash. Factorization is independent of the active-space selection, so it can run after step 02 if desired. Estimate memory and disk before launching; the allocation holds the packed AO factors in memory.

Scripts: [ao_cholesky.py](ao_cholesky.py), [check_ao_factors.py](check_ao_factors.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 07_ao_cholesky/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 07_ao_cholesky/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 16 CPU(s), 64 GiB RAM, walltime 04:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `geometry.xyz`
- `parent/scf.chk`
- `inputs/basis.json`

Key generated outputs (other diagnostic/state files are also written):

- `ao_factorization/metadata.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

Tight rank 9218; loose rank 7545. AO packed factors occupy approximately 15.07 GiB and are omitted from Git. Metadata and sampled checks are retained. A dense retained-space ERI tensor alone would need roughly 1.06 TiB (618^4 float64 elements), so the factor representation is essential.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/ao_factorization/metadata.json](results/ao_factorization/metadata.json)
- [results/validation/ao_factor_samples.json](results/validation/ao_factor_samples.json)

## Reuse for a similar task

Adjust rank capacity and memory estimates for the new AO basis; validate both residual and downstream observable errors. The factor threshold is not the final many-body error tolerance.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
