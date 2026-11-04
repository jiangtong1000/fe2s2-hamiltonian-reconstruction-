# Provenance

- `original_basis_and_method.txt`: the implemented method and settings of the completed 2026-10-05/06 calculation (geometry source, parent SCF, relativity, localization, frozen core, DMRG and factorization choices). Paths in it refer to the original calculation directory, which is not part of this repository.
- `repository_checksums.json`: SHA-256 and size of every committed file. `python tools/check_repository.py` verifies it; `--update-inventory` regenerates it after a reviewed edit.

Recorded JSON under `*/results/` keeps the metadata of the original runs, including relative paths in the original directory and Slurm job IDs. These identify the completed calculation; they are not runnable paths here. The historical reference FCIDUMP is not redistributed; step 01 fetches it from its pinned upstream commit and verifies the recorded SHA-256.
