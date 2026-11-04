"""Hash durable outputs; report excluded rebuildable scratch without deleting it."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_shared"))
from workflow_paths import DATA_ROOT, CODE_ROOT
import datetime
import hashlib
import importlib.metadata
import json
import subprocess


def sha(path: Path) -> str:
    """Hash even multi-gigabyte factors without loading them into memory."""
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    required = [
        "parent/density_alpha_ao.npy",
        "parent/density_beta_ao.npy",
        "parent/mo_coeff.npy",
        "parent/parent_natural_occupations.npy",
        "cas20_provisional/rdm1.npy",
        "orbitals/U_natural_20.npy",
        "orbitals/natural_orbital_diagnostics.json",
        "full_hamiltonian/h1.npy",
        "full_hamiltonian/chol.npy",
        "full_hamiltonian/h0.txt",
        "full_hamiltonian/index_maps.json",
        "orbitals/trial_embeddings.json",
        "cas20/rdm1.npy",
        "cas20/rdm2.npy",
        "cas14/rdm1.npy",
        "cas14/rdm2.npy",
        "validation/validation_metrics.json",
        "validation/cas14_fci_comparison.json",
        "validation/scientific_assessment.json",
    ]
    for name in required:
        if not (DATA_ROOT / name).is_file():
            raise FileNotFoundError(name)
    for name in ["reference_zhendong", "cas20_provisional", "cas20", "cas14"]:
        if not any((DATA_ROOT / name / "native_state").glob("*GS*")):
            raise FileNotFoundError(f"{name}/native_state: no GS state files")
    validation = json.loads(
        (DATA_ROOT / "validation/validation_metrics.json").read_text()
    )
    if validation["status"] != "numerical_checks_passed":
        raise RuntimeError("Numerical checks have not passed")
    fci = json.loads((DATA_ROOT / "validation/cas14_fci_comparison.json").read_text())
    if fci["status"] != "passed":
        raise RuntimeError("Independent FCI comparison has not passed")
    records = {}
    excluded = []
    for path in sorted(DATA_ROOT.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(DATA_ROOT)
        name = str(rel)
        # Checkpoint tensors and FCI vectors remain durable; omit intermediate MPO/operator scratch.
        scratch = any(
            "scratch" in part
            or part.startswith(("native_state_before", "before_refinement_stage"))
            for part in rel.parts
        )
        scratch |= "native_state" in rel.parts and "GS" not in path.name
        if name == "run_manifest.json" or scratch:
            excluded.append({"path": name, "size_bytes": path.stat().st_size})
            continue
        records[name] = {"sha256": sha(path), "size_bytes": path.stat().st_size}
    packages = {}
    for name in [
        "numpy",
        "scipy",
        "pyscf",
        "block2-mpi",
        "h5py",
        "matplotlib",
        "PyYAML",
    ]:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    revision = subprocess.run(
        ["git", "-C", str(CODE_ROOT), "rev-parse", "HEAD"],
        text=True,
        capture_output=True,
    )
    result = {
        "schema_version": 1,
        "status": "archived_with_scientific_review_pending",
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_revision": revision.stdout.strip() if revision.returncode == 0 else None,
        "code_sha256": {
            str(p.relative_to(CODE_ROOT)): sha(p)
            for p in sorted(CODE_ROOT.rglob("*.py"))
            if "work" not in p.relative_to(CODE_ROOT).parts
            and "provenance" not in p.relative_to(CODE_ROOT).parts
        },
        "python": sys.version,
        "packages": packages,
        "parent_natural_orbitals": {
            "source": "spin-summed high-spin UKS 1-RDM",
            "density_artifacts": [
                "parent/density_alpha_ao.npy",
                "parent/density_beta_ao.npy",
            ],
            "coefficients": "parent/mo_coeff.npy",
            "occupations": "parent/parent_natural_occupations.npy",
            "run_provenance": "run_logs",
        },
        "cas20_natural_orbitals": {
            "source": "singlet CAS(30e,20o) DMRG 1-RDM",
            "rotation_source_stage": "provisional_direct_AO",
            "density": "cas20_provisional/rdm1.npy",
            "rotation": "orbitals/U_natural_20.npy",
            "occupations": "orbitals/natural_orbital_diagnostics.json",
            "production_density": "cas20/rdm1.npy",
            "run_provenance": "run_logs",
        },
        "partition": "full_hamiltonian/index_maps.json",
        "trial_embeddings": "orbitals/trial_embeddings.json",
        "trial_fixed_occupations_restrict_full_propagation": False,
        "artifacts": records,
        "excluded_rebuildable_or_self_files": excluded,
        "note": "All listed hashes were computed in this run. Review scientific_assessment.json; archive creation does not accept CAS14 for an application. No files were deleted.",
    }
    (DATA_ROOT / "run_manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"Archived {len(records)} files; {len(excluded)} scratch/self files excluded; no deletion."
    )


if __name__ == "__main__":
    main()
