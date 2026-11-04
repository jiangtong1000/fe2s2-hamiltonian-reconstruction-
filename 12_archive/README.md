# 12_archive: Reproducibility manifest and storage policy

Record a verifiable archive while keeping Git compact and heavy state data separate.

## Method and scripts

`archive_run.py` hashes durable artifacts in a new work directory, records package versions and code hashes, and writes `run_manifest.json`. It preserves the two distinct natural-orbital provenance sections and the semantic partition/trial maps. It requires numerical validation and current native GS states, but does not automatically approve compression for an application.

It excludes regenerable operator scratch and historical checkpoint backups from its durable inventory; it deletes nothing. Run this step after inspecting the scientific assessment. The portable manifest is JSON rather than the original hand-curated YAML, and does not invent historical job IDs or transplant the archived completion status.

Scripts: [archive_run.py](archive_run.py). `run.py` is the step entry point; `step.json` records commands, prerequisites, and suggested resources.

## Run

From the repository root, inspect the command first:

```bash
python 12_archive/run.py --workdir /scratch/fe2s2_run --dry-run
```

Then run inside an appropriate compute allocation:

```bash
python 12_archive/run.py --workdir /scratch/fe2s2_run
```

Actions:

- `default`.

Use `--action NAME` to select an alternative. For refinement actions, `--bonds ...` overrides the saved defaults. Suggested resources: 2 CPU(s), 8 GiB RAM, walltime 02:00:00. These are planning estimates, not performance guarantees. See [environment instructions](../environment/README.md) for batch execution.

Prerequisites (paths relative to the work directory):

- `validation/validation_metrics.json`
- `validation/compression_assessment.json`

Key generated outputs (other diagnostic/state files are also written):

- `run_manifest.json`

New calculations write to the work directory, never to the committed `results/` folder. Run step 01 once to initialize a new directory. Original metadata inside `results/` may refer to paths in the original archive; those are provenance, not runnable paths in this repository.

## Recorded result

The original manifest and completion audits are retained under `results/` and `../provenance/`. The original audit covered 553 hashed files totaling 44,999,098,931 bytes, plus its own report entry. `../provenance/omitted_artifacts.json` identifies large omitted files, original checksums, and rebuild steps. The original local calculation archive remains intact.

These are results of the completed 2026-10-05/06 calculation. Packaging checks do not constitute a new molecular rerun.

Retained evidence:

- [results/validation/artifact_audit.json](results/validation/artifact_audit.json)
- [results/validation/completion_audit.json](results/validation/completion_audit.json)

## Reuse for a similar task

Keep immutable inputs, small diagnostic results, code, and provenance in Git. Retain expensive MPS/FCI checkpoints in durable external storage; factors and operator scratch are better rebuild candidates. A compact Git clone alone does not contain the full VCE Hamiltonian and state data.

[Workflow overview](../README.md) · [Storage policy](../STORAGE.md)
