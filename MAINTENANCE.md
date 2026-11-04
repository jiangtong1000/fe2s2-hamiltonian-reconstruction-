# Maintaining this workflow

Keep the twelve numbered folders: they are the scientific narrative and public entry points. A step stays understandable through its README, script, `step.json`, and recorded results. Prefer a small change in the existing file over a new abstraction or dependency.

## Where to change code

| Concern | File |
|---|---|
| Commands, prerequisites, refinement choices, resource suggestions | Each step's `step.json` |
| CLI planning, work-directory protection, sequential execution and logs | `_shared/step_runner.py` |
| Repository and calculation roots | `_shared/workflow_paths.py` |
| Physical fixed-core and embedded-energy contractions | `_shared/integrals.py` |
| Shared singlet solver and state/RDM conventions | `_shared/run_dmrg.py` |
| Molecule-specific chemistry and numerical validation | The relevant numbered scientific script |
| Source/data integrity | `tools/check_repository.py` and provenance inventories |

No framework or package installation is required for step discovery, dry runs, or packaging tests. The small per-step `run.py` files are deliberate: they keep every numbered folder directly executable. The scientific scripts share only code they actually use.

## Code and provenance rules

- Use descriptive names, ordinary control flow, module/function documentation where conventions matter, and Black formatting. Avoid semicolon-packed statements.
- Keep the scientific assertions, directory protections, failed-run logs, and independent validation routes. Do not run with Python `-O`, which disables assertions.
- Treat `*/results/`, `01_inputs/inputs/`, `provenance/original_scripts/`, and `provenance/submitted/` as immutable historical evidence. The formatter excludes results and provenance. Do not rewrite old results to make a changed algorithm appear validated.
- Keep formula/order changes separate from formatting changes. Any scientific change needs the relevant numerical validation and fresh output provenance.
- Preserve parent-UKS versus singlet-CAS natural-orbital provenance, the block2 2-RDM index convention, the scalar energy definition, and the full F/I/K6/A14/V semantics.
- Preserve original source hashes. `portable_script_changes.json` records original and current portable hashes; add an accurate change description when editing those scripts. The current full-file inventory is separate from the historical calculation manifest.

## Local checks

From the repository root, using Python 3.11 and an available Black installation:

```bash
# One formatter worker avoids unnecessary processes on shared/login nodes.
python -m black --workers 1 .
python -m black --workers 1 --check .
python -m unittest discover -s tests -v
```

After inspecting the diff, deliberately refresh the current repository inventory:

```bash
python tools/check_repository.py --update-inventory
python tools/check_repository.py
git diff --check
git status --short
```

The inventory excludes itself, Git metadata, caches, and generated `work/` data. It still verifies the original-copy and input hashes before updating, so regenerating the inventory cannot silently bless edits to those protected files. Record new validation results separately from `provenance/packaging_validation.json`, which describes the initial packaging.

For numerical changes, use the installed chemistry environment and the appropriate compute allocation. Step 10's `smoke` action supplies small contraction/factorization/projector checks; complete molecular checks belong in Slurm. Packaging checks alone never certify a new molecular result.

## Code style used here

Reuse existing code and standard-library facilities, avoid speculative infrastructure, and retain correctness safeguards. During the maintenance cleanup, scientific operations were checked for syntax-tree equivalence apart from imports and documentation; the runner and streaming hash helper received explicit implementation changes.
