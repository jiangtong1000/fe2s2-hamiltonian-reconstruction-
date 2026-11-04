# 08_full_hamiltonian: Full retained Hamiltonian and both CAS Hamiltonians

Build a backend-neutral full Hamiltonian and derive both final CAS Hamiltonians from it.

## Method and scripts

`build_production.py` transforms AO Cholesky factors to all 618 retained orbitals and contracts **only F44**. It stores h1, the nuclear-plus-frozen-core scalar h0, and factors L[g,p,q], with chemists’ ERIs (pq|rs)=Σg L[g,p,q]L[g,r,s]. The operator is physical: no DFT XC potential or AFQMC-specific shift is included.

Both final CAS files come from this same serialized representation: freeze I to obtain CAS(30e,20o), or I+K6 to obtain CAS(18e,14o). The loose-rank representation includes its own consistently evaluated frozen-core contraction. Full cross-block integrals, explicit F/I/K6/A14/V roles, electron counts, and trial maps survive.

Scripts: [build_production.py](build_production.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 08_full_hamiltonian/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 08_full_hamiltonian/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 16 CPU(s), 64 GiB RAM, walltime 04:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `orbitals/index_maps.json`
- `ao_factorization/chol_ao_packed.npy`

Key generated outputs (other diagnostic/state files are also written):

- `full_hamiltonian/schema.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

Full factors have shape (9218,618,618), approximately 26.23 GiB; the file is omitted from Git. Both production FCIDUMPs, h0, electron count, and schema are retained. CAS20’s derived scalar is −4976.396220614395 Eh; it is not the historical scalar.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/cas14/FCIDUMP](results/cas14/FCIDUMP)
- [results/cas20/FCIDUMP](results/cas20/FCIDUMP)
- [results/full_hamiltonian/h0.txt](results/full_hamiltonian/h0.txt)
- [results/full_hamiltonian/nelec.txt](results/full_hamiltonian/nelec.txt)
- [results/full_hamiltonian/schema.json](results/full_hamiltonian/schema.json)

## Reuse for a similar task

Do not flatten the semantic partition or eliminate I/K6/V based on trial occupations. Future VCE/active-space AFQMC needs these blocks even though backend code is outside this repository’s completed phase.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
