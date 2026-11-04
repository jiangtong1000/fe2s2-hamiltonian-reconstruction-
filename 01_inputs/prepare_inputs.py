"""Install the committed inputs and the external reference FCIDUMP; verify SHA-256s.

The historical FCIDUMP is third-party data and is not part of this repository. It is
taken from the local file named by FE2S2_REFERENCE_FCIDUMP when that variable is set,
and otherwise downloaded from the pinned upstream commit in inputs/external_inputs.json.
Either way it is refused unless its SHA-256 matches; inputs are never substituted.
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_shared"))
from workflow_paths import DATA_ROOT, CODE_ROOT
import hashlib
import json
import os
import shutil
import urllib.request

INPUTS = CODE_ROOT / "01_inputs/inputs"


def reference_fcidump(dest):
    """Copy or download the reference FCIDUMP to dest and verify its SHA-256."""
    external = json.loads((INPUTS / "external_inputs.json").read_text())
    spec = external["reference.FCIDUMP"]
    local = os.environ.get(spec["local_copy_variable"])
    partial = dest.with_name(dest.name + ".partial")
    if local:
        shutil.copyfile(Path(local).expanduser(), partial)
        origin = str(Path(local).expanduser().resolve())
    else:
        with urllib.request.urlopen(spec["url"], timeout=60) as response:
            partial.write_bytes(response.read())
        origin = spec["url"]
    digest = hashlib.sha256(partial.read_bytes()).hexdigest()
    if digest != spec["sha256"]:
        partial.unlink()
        raise ValueError(
            f"Reference FCIDUMP SHA-256 mismatch from {origin}: "
            f"{digest} != {spec['sha256']}"
        )
    partial.replace(dest)
    return origin, digest


def main():
    inventory = json.loads((INPUTS / "checksums.json").read_text())
    for name, digest in inventory.items():
        assert hashlib.sha256((INPUTS / name).read_bytes()).hexdigest() == digest, name
    mapping = {
        "geometry.xyz": "geometry.xyz",
        "basis.json": "inputs/basis.json",
        "geometry_provenance.json": "sources/geometry_provenance.json",
        "reference_provenance.json": "reference_zhendong/provenance.json",
        "reference_initial_occupations.npy": "reference_zhendong/initial_occupations.npy",
    }
    fcidump = DATA_ROOT / "reference_zhendong/FCIDUMP"
    for dest in [fcidump, *(DATA_ROOT / target for target in mapping.values())]:
        if dest.exists():
            raise FileExistsError(dest)
    fcidump.parent.mkdir(parents=True, exist_ok=True)
    origin, digest = reference_fcidump(fcidump)
    for source, target in mapping.items():
        dest = DATA_ROOT / target
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(INPUTS / source, dest)
    record = dict(inventory, **{"reference.FCIDUMP": digest})
    (DATA_ROOT / "inputs/input_checksums.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n"
    )
    (DATA_ROOT / "reference_zhendong/FCIDUMP_source.json").write_text(
        json.dumps({"origin": origin, "sha256": digest}, indent=2) + "\n"
    )
    print(
        "Verified and installed the geometry, explicit basis, and historical FCIDUMP "
        f"(from {origin})."
    )


if __name__ == "__main__":
    main()
