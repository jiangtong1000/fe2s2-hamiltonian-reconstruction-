"""Verify compact repository payload and runnable step metadata (stdlib only)."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "provenance/repository_checksums.json"


def payload():
    for p in sorted(ROOT.rglob("*")):
        rel = p.relative_to(ROOT)
        if any(
            x in {".git", "work", "__pycache__", ".pytest_cache"} for x in rel.parts
        ):
            continue
        if p.is_file() and p != INVENTORY:
            yield p


def records():
    return {
        str(p.relative_to(ROOT)): {
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            "size_bytes": p.stat().st_size,
        }
        for p in payload()
    }


def verify_structure():
    steps = sorted(p for p in ROOT.iterdir() if p.is_dir() and p.name[:2].isdigit())
    assert len(steps) == 12
    for i, step in enumerate(steps, 1):
        assert step.name.startswith(f"{i:02d}_")
        assert (step / "README.md").is_file() and (step / "run.py").is_file()
        spec = json.loads((step / "step.json").read_text())
        if spec.get("recorded_results", True):
            assert any((step / "results").rglob("*")), step.name
        assert spec["step"] == step.name and "default" in spec["actions"]
        for action in spec["actions"].values():
            for command in action["commands"]:
                assert (ROOT / command["script"]).is_file(), command
                assert not Path(command["cwd"]).is_absolute(), command
    for p in payload():
        if p.suffix == ".py":
            ast.parse(p.read_text(), filename=str(p))
        if p.suffix == ".json":
            json.loads(p.read_text())
        assert p.stat().st_size < 10 * 1024**2, f"Unexpected large Git file: {p}"
        if p.suffix == ".md":
            for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", p.read_text()):
                if "://" in target or target.startswith("#"):
                    continue
                target = target.split("#")[0]
                assert (
                    p.parent / target
                ).exists(), f"Broken local link in {p}: {target}"
    inputs = ROOT / "01_inputs/inputs"
    for name, digest in json.loads((inputs / "checksums.json").read_text()).items():
        assert hashlib.sha256((inputs / name).read_bytes()).hexdigest() == digest, name
    return len(steps)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--update-inventory",
        action="store_true",
        help="Intentionally regenerate hashes after reviewing repository edits",
    )
    args = parser.parse_args()
    count = verify_structure()
    current = records()
    if args.update_inventory:
        INVENTORY.write_text(json.dumps(current, indent=2, sort_keys=True) + "\n")
    else:
        expected = json.loads(INVENTORY.read_text())
        assert (
            expected == current
        ), "Payload differs from inventory; inspect changes before updating checksums"
    print(
        f'PASS: {count} steps; {len(current)} checked files; {sum(x["size_bytes"] for x in current.values())/1024**2:.2f} MiB payload; syntax, links, input hashes, and action graph valid.'
    )


if __name__ == "__main__":
    main()
