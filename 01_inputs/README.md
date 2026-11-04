# 01_inputs: Inputs and provenance

Establish the geometry, explicit AO basis, electron count, and historical comparison before doing any calculation.

## Method and scripts

`prepare_inputs.py` verifies SHA-256 checksums and copies the committed inputs into a new work directory. `inputs/geometry.xyz` is the oxidized geometry from Sharma et al. 2014 SI Table 1. `inputs/basis.json` is the explicit TZP-DKH definition recorded by the completed parent calculation; it avoids silently accepting a changed basis library.

The historical FCIDUMP is third-party data and is not redistributed. `prepare_inputs.py` downloads it from the pinned upstream commit recorded in `inputs/external_inputs.json` and refuses it unless its SHA-256 matches. On a machine without internet access, download it once and set `FE2S2_REFERENCE_FCIDUMP=/path/to/fe2s2`; the same check applies. The historical FCIDUMP has a zero scalar. Its external core energy is −4976.26532397 Eh, supplied separately by the pinned upstream README. Add it once when comparing total energies; never copy it into our reconstructed Hamiltonian.

Scripts: [prepare_inputs.py](prepare_inputs.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 01_inputs/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 01_inputs/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 1 CPU(s), 2 GiB RAM, walltime 00:10:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- No generated prerequisites; committed offline inputs are sufficient.

Key generated outputs (other diagnostic/state files are also written):

- `geometry.xyz`
- `reference_zhendong/FCIDUMP`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

24 atoms; 186 electrons; charge −2. Source coordinates have 0.01 Å precision. Fe–Fe = 2.68648469 Å. The historical FCIDUMP SHA-256 is `95d8786af06eeea2107e19ffd98c66a6ca97fc8c9864175a4f6d64512b6f2df9`.

The source PDF is omitted; its URL, table/page, and checksum are in `inputs/geometry_provenance.json`.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/build/resource_estimates.json](results/build/resource_estimates.json)
- [results/sources/geometry_provenance.json](results/sources/geometry_provenance.json)

## Reuse for a similar task

Change the molecule, geometry provenance, charge/spin, explicit basis, and reference attribution together. Downstream scripts contain Fe2S2-specific assumptions; replacing the XYZ alone is insufficient.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
