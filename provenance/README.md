# Provenance and portability changes

`original_manifest.yaml` preserves the completed calculation's settings, artifacts and job records. It is a historical record, not an instruction to restart cluster jobs or a claim that a new run has passed. For the public release, local filesystem paths, internal review notes and entries outside the scope of this repository were removed from seven retained records; `retained_files.json` lists their original and current hashes and what changed. The historical reference FCIDUMP is not redistributed; step 01 fetches it from its pinned upstream commit and verifies the recorded SHA-256. Original paths are relative to the original `fe2s2_full/` archive unless explicitly absolute.

`original_scripts/` contains the unmodified scientific source versions used as the packaging baseline. `submitted/` preserves available submitted snapshots, including earlier solver/restart history. These historical scripts can contain hard-coded paths and job IDs; use the numbered portable step entry points for new work.

`portable_script_changes.json` records provenance for relocated scientific scripts. Their formulas and contraction conventions are preserved. Edits make data roots configurable, resolve shared imports, use committed explicit basis definitions, route helper subprocesses to their new locations, and make incomplete/failed full validation return a failure exit status. Packaging adds new step runners, offline input preparation, an assessment summary, and a generic archive writer. These new helpers do not inherit the original molecular validation merely by being copied into the repository.

The corrected shared DMRG driver resets to two-site optimization at each bond stage. Its fresh-run schedule cannot recreate every historical restart, and the recorded original results are not presented as newly generated outputs of the wrappers.

`retained_files.json` gives hashes for direct copies. `omitted_artifacts.json` records original manifest hashes and rebuild routes for large files not directly copied. RDM and orbital NPZ bundles are lossless reorganizations of original small arrays; the repository checker verifies their recorded checksums, and packaging tests compared them with the originals. `repository_checksums.json` covers the current repository payload (excluding itself, ignored work, Git metadata, and Python caches).

Run `python tools/check_repository.py` to verify the payload and command graph. Packaging tests cover syntax, input checksums, output-directory protection, all actions' dry runs, and documented file links. They do not rerun SCF, DMRG, factorization, or molecular integral checks.

[Packaging validation](packaging_validation.json) records the checks performed during this reorganization. [Compressed array sources](compressed_array_sources.json) maps each NPZ key to its original array and records exact equality.

[Maintenance cleanup validation](cleanup_validation.json) records the later maintenance cleanup, its computational-equivalence checks, and runner regression tests. Historical reports above remain unchanged.
