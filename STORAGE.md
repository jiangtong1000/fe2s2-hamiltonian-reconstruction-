# Git contents and external calculation data

The compact repository preserves small inputs, scripts, final CAS FCIDUMPs, compressed 1-/2-RDM arrays, the final orbital coefficient/rotation bundle, maps, metadata, convergence histories, validation summaries, and PDF figures. `provenance/retained_files.json` maps exact copied files to original paths and hashes; compressed bundles are described below. `provenance/repository_checksums.json` checks the entire committed payload apart from itself and Git metadata.

The original calculation directory is not part of this repository. The original manifest records full artifact hashes and original locations; `provenance/omitted_artifacts.json` lists large omitted artifacts with rebuild steps. These pointers document provenance, not guaranteed future availability.

| Omitted class | Why omitted from Git | Rebuild / retention |
|---|---|---|
| AO packed Cholesky, ~15.07 GiB | Large deterministic intermediate | Step 07, using exact geometry/basis |
| Full retained factors, ~26.23 GiB | Large deterministic Hamiltonian representation | Steps 06–08; keep exact coefficients/maps if regenerating from the archived basis |
| Parent checkpoint and most AO/MO arrays | Reproducible but machine/version-sensitive numerical result | Step 02; preserve externally when exact historical coefficients matter |
| Native MPS states and operator scratch | States can be expensive to obtain; scratch is rebuildable | Keep GS checkpoint/state metadata externally; regenerate MPO scratch. Rerunning steps 04/05/09 is possible but costly |
| CAS14 FCI vector | Large solver state | Step 10; preserve externally if repeated FCI-state analyses are expected |
| Refinement backups and verbose scheduler logs | Large/redundant history | Compact convergence reports and submitted-source snapshots retained; original files remain local |
| Source PDFs | Readily obtainable from recorded sources | Source URLs, identifying table/page, and hashes retained |
| Historical reference FCIDUMP | Third-party data, not redistributed | Step 01 fetches it from the pinned upstream commit and verifies its SHA-256 |

“Not committed” does not mean “deleted.” In particular, expensive many-body states should not be treated as disposable just because they are reproducible. Exact historical recomputation may depend on the original runtime and numerical path.

## Compressed arrays

Each `*_rdms.npz` in steps 04/05/09 uses original `.npy` stems as keys (for example `rdm1`, `rdm2`, and `rdm1_stage5`). Load with `numpy.load(..., allow_pickle=False)`. Compression is lossless; it does not truncate states or occupations. No claim is made that an RDM bundle replaces an MPS for arbitrary higher-order observables or overlaps.

`06_natural_orbitals/results/final_basis.npz` retains the final full coefficients, frozen-core coefficients, initial CAS20 coefficients, and the natural-orbital rotation. Inspect `.files` for exact keys. The orbital index maps and source density provenance are stored alongside it. This small bundle preserves the chosen gauge even though the parent checkpoint is omitted.

## Full-Hamiltonian rebuild

The runnable route is steps 01–08 in a fresh work directory, followed by state and validation steps. For a new parent calculation, review chemical selection anew and expect numerically different localized coefficients. To reproduce the archived representation exactly, retain the original parent physical hcore/overlap/metadata and orbital arrays externally, or restore them from the original archive before using the scientific builder. The compact clone alone is not a byte-identical restore of that full archive.

A full VCE/active-space AFQMC handoff additionally needs the factor arrays, semantic maps, trial state data, and conventions. Those remain in the original calculation archive or can be regenerated and revalidated. The compact Git repository intentionally does not claim to be that full handoff.

## New runs

Put generated data under ignored `work/` or external scratch. Step 12 writes a fresh hash manifest; it does not delete files or add data to Git. Commit reviewed small summaries deliberately, with run provenance. Check `git status` and file sizes before every commit. No Git LFS or remote storage service is required for this reference repository.
