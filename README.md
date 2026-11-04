# Fe2S2 Hamiltonian reconstruction: a reusable worked example

A step-by-step record of reconstructing the oxidized [Fe2S2(SCH3)4]2− Hamiltonian, selecting CAS(30e,20o), and assessing a nested CAS(18e,14o). Each numbered folder contains a runnable entry point, implementation scripts or an explicit shared-script link, a README, and compact recorded results. Original calculation: 2026-10-05/06.

This repository preserves the **method and evidence**, including the finding that six-orbital freezing costs about 34.17 mEh. Energy agreement did not drive orbital selection. The reconstructed Hamiltonian is chemically consistent with the historical model. With the reference orbitals published on 2026-10-10, it reproduces the historical core energy and active one-electron integrals, and the reconstructed CAS(30e,20o) spans essentially the same space with a 1.43 mEh lower CASCI energy; see [docs/active_space_comparison](docs/active_space_comparison/README.md).

## Steps

| Step | Purpose | Main recorded outcome |
|---|---|---|
| [01_inputs](01_inputs/README.md) | Geometry, basis, reference provenance | Sharma 2014 oxidized geometry; explicit TZP-DKH basis |
| [02_parent_uks](02_parent_uks/README.md) | High-spin BP86/sf-X2C parent and natural orbitals | 186 electrons, multiplicity 11, 662 spherical AOs |
| [03_localize_select](03_localize_select/README.md) | Split PM localization and chemical selection | 10 Fe 3d + 6 bridge S 3p + 4 terminal S 3p |
| [04_reference_dmrg](04_reference_dmrg/README.md) | Historical singlet DMRG/RDM benchmark | Eactive = −116.605608594803 Eh |
| [05_provisional_cas20](05_provisional_cas20/README.md) | Direct-AO CAS20 and provisional correlated density | Density for the fixed compression basis |
| [06_natural_orbitals](06_natural_orbitals/README.md) | Correlated NOs, nested spaces, trial maps | F44 separate; retained I34 / K6 / A14 / V564 |
| [07_ao_cholesky](07_ao_cholesky/README.md) | Factor the AO Coulomb operator | Rank 9218 at residual threshold 1e-9 Eh |
| [08_full_hamiltonian](08_full_hamiltonian/README.md) | Full Hamiltonian and both production CAS files | 618 orbitals, 98 electrons; all cross-block couplings |
| [09_production_states](09_production_states/README.md) | Final CAS20 and CAS14 singlets | CAS20 total = −5092.872362559896 Eh |
| [10_independent_validation](10_independent_validation/README.md) | Integral, embedding, virtual-block, FCI checks | Numerical checks passed; CAS14 FCI agrees within ~8.7e-10 Eh |
| [11_assessment](11_assessment/README.md) | Historical occupations and compression cost | ΔEfreeze = 34.1692 mEh; K6-double weight = 0.940363 |
| [12_archive](12_archive/README.md) | Artifact inventory and reproducibility | Checksums, provenance, external-data/rebuild policy |

Steps 04 and 07 can run independently once their inputs exist. The numbered order is a convenient serial route. Step 06 requires a sufficiently converged provisional density; step 11 requires final state and validation results. Review diagnostics between steps rather than blindly chaining jobs.

## Inspect or rerun

Read the numbered READMEs and their committed `results/` without running any chemistry. No external archive is needed to inspect the reconstructed FCIDUMPs, RDMs, final basis, plots, or recorded checks.

```bash
# Packaging verification: standard-library checks; no molecular computation.
python tools/check_repository.py
python -m unittest discover -s tests -v

# Initialize a NEW directory: copy checksum-verified inputs, fetch the reference FCIDUMP.
python 01_inputs/run.py --workdir /scratch/fe2s2_example

# Preview a step without importing chemistry packages or creating files.
python 02_parent_uks/run.py --workdir /scratch/fe2s2_example --dry-run

# On a configured compute node:
python 02_parent_uks/run.py --workdir /scratch/fe2s2_example
```

Use the same `--workdir` for subsequent steps, or export `FE2S2_WORKDIR`. The default is this repository's ignored `work/` directory. An existing unmarked calculation directory is refused, so earlier calculations are never overwritten. See [environment/README.md](environment/README.md) for dependencies and Slurm commands. Do not run the molecular calculations on a login node.

DMRG schedules are starting points. The completed production CAS20 needed refinement through bond dimension 5000. A new run must demonstrate its own energy, occupation, and subspace convergence; bitwise reproduction is not promised across machines or numerical libraries.

## What is preserved

Small inputs and results are committed. Compressed NPZ files preserve the saved RDM arrays and final coefficient/rotation matrices. The production and historical FCIDUMPs, chemical reviews, index maps, convergence histories, figures, and independent validation metrics remain directly inspectable. [STORAGE.md](STORAGE.md) explains omitted arrays, reconstruction routes, and which expensive states deserve external retention.

`results/` always means **recorded original evidence**, not output from the portable wrappers. Copied JSON/YAML may contain original relative paths and Slurm job IDs. Those identify the completed calculation; new runs record their own paths, command hashes, logs, versions, and states. [provenance/README.md](provenance/README.md) explains source copies and the limited portability edits. Original calculation files were left intact.

Two natural-orbital constructions remain explicitly separate:

```yaml
parent_natural_orbitals:
  source: spin-summed high-spin UKS 1-RDM
cas20_natural_orbitals:
  source: singlet CAS(30e,20o) DMRG 1-RDM
```

The CAS20 rotation uses the provisional direct-AO density. Final production RDMs validate that fixed basis and assess compression. Only genuine F44 core orbitals are contracted out of the full Hamiltonian; trial-fixed I/K6 and trial-empty V do not restrict future propagation. AFQMC adapters, propagation, and truncated MSD export were outside this phase.

For code changes and routine checks, see [MAINTENANCE.md](MAINTENANCE.md).

## Limits and adapting the example

This is a worked Fe2S2 workflow, not a general-purpose automatic active-space selector. Atom indices, electron counts, core shells, chemical assignment rules, orbital block sizes, reference core energy, and some validation assertions are molecule-specific. For a new task, change and review them together. Keep an independent direct-integral route and explicit external-virtual checks.

Matching total energies or natural occupations does not prove historical Hamiltonian identity. Our independent integral routes share the PySCF AO engine and input model. CAS14’s sector weight is not fidelity with its reoptimized state and does not predict future AFQMC bias. See the individual steps for the actual evidence and limitations.

## Sources and attribution

- [Sharma et al. 2014 SI](https://authors.library.caltech.edu/records/xzxp2-7q641/files/nchem.2041-s1.pdf?download=1): oxidized geometry, Table 1. Source hash and extraction details are committed.
- [Li–Chan 2017](https://arxiv.org/pdf/1703.04789): later natural-orbital/split-localization route.
- [Lee et al. 2023 SI](https://rcastoragev2.blob.core.windows.net/0a8f19fd29e3811703d0eafa850f0e3a/41467_2023_37587_MOESM1_ESM.pdf): model construction and references.
- [Zhendong Li's pinned reference README](https://github.com/zhendongli2008/Active-space-model-for-Iron-Sulfur-Clusters/blob/84433831a680f70665cf43fc8692e0210da92be3/README.md): published active/core energies; its FCIDUMP is not redistributed here; step 01 fetches it from the pinned commit and verifies its SHA-256.
- Sharma–Chan 2012, DOI [10.1063/1.3695642](https://doi.org/10.1063/1.3695642): chemical 20-orbital composition/order. Its older SVP model is not the later TZP-DKH Hamiltonian.

## License

The code and the data generated by this workflow are released under the [MIT License](LICENSE). Third-party data and source documents keep their original attribution and terms: the historical reference FCIDUMP is not included and is fetched from its upstream repository, and the cited publications are linked, not redistributed.
