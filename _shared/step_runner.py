"""Plan and run numbered steps without importing scientific dependencies.

Only execution creates files. Work-directory checks protect archived results;
command logs record failed runs as well as successful ones.
"""

import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from workflow_paths import CODE_ROOT as REPO

MARKER = ".fe2s2_workdir.json"
WORKFLOW_NAME = "fe2s2_reference"


def build_commands(action: dict, work: Path, bonds: list[int]) -> list[tuple]:
    """Resolve the selected action into argv, working directory, and source path."""
    commands = []
    for command in action["commands"]:
        script = REPO / command["script"]
        argv = [sys.executable, str(script), *command.get("args", [])]
        if command.get("refine"):
            argv.extend(["--refine-bonds", *map(str, bonds)])
        commands.append((argv, work / command.get("cwd", "."), script))
    return commands


def prepare_workdir(work: Path, step: Path, spec: dict, action: dict) -> None:
    """Initialize a fresh directory or validate a previously initialized run."""
    if work == REPO or (
        work.is_relative_to(REPO) and not work.is_relative_to(REPO / "work")
    ):
        raise ValueError(
            "Use repo/work or a separate scratch directory, not a checked-in step directory"
        )
    marker_path = work / MARKER
    if not marker_path.exists():
        if step.name != "01_inputs":
            raise ValueError("Initialize the work directory with step 01 first")
        if work.exists() and any(work.iterdir()):
            raise ValueError("Refusing a nonempty unmarked work directory")
        work.mkdir(parents=True, exist_ok=True)
        marker_path.write_text(
            json.dumps(
                {
                    "workflow": WORKFLOW_NAME,
                    "created_utc": datetime.datetime.now(
                        datetime.timezone.utc
                    ).isoformat(),
                },
                indent=2,
            )
            + "\n"
        )
    if json.loads(marker_path.read_text()).get("workflow") != WORKFLOW_NAME:
        raise ValueError("Unrecognized work-directory marker")

    missing = [
        name
        for name in action.get("requires", spec["requires"])
        if not (work / name).exists()
    ]
    if missing:
        raise ValueError("Missing prerequisite(s): " + ", ".join(missing))
    existing = [name for name in action.get("protect", []) if (work / name).exists()]
    if existing:
        raise ValueError(
            "Outputs already exist; use a new work directory or the documented "
            "continuation action: " + ", ".join(existing)
        )
    for name in [
        "inputs",
        "sources",
        "parent",
        "orbitals",
        "reference_zhendong",
        "cas20_provisional",
        "ao_factorization",
        "full_hamiltonian",
        "cas20",
        "cas14",
        "cas20_loose",
        "cas14_loose",
        "validation",
        "run_logs",
    ]:
        (work / name).mkdir(exist_ok=True)


def execute_commands(
    step: Path, action_name: str, commands: list[tuple], work: Path, env: dict[str, str]
) -> None:
    """Run commands in order and save provenance even when a command fails."""
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    record = {
        "step": step.name,
        "action": action_name,
        "workdir": str(work),
        "started_utc": stamp,
        "slurm_job_id": env.get("SLURM_JOB_ID"),
        "commands": [],
        "environment": {
            key: env.get(key)
            for key in [
                "OMP_NUM_THREADS",
                "DMRG_STACK_GB",
                "VALIDATION_MEMORY_MB",
                "FE2S2_RUNTIME",
            ]
        },
        "status": "failed",
    }
    record_path = work / "run_logs" / f"{stamp}_{step.name}_{action_name}.json"
    try:
        for number, (argv, cwd, script) in enumerate(commands):
            log = work / "run_logs" / f"{stamp}_{step.name}_{number}.log"
            item = {
                "argv": argv,
                "cwd": str(cwd),
                "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
                "log": str(log.relative_to(work)),
            }
            record["commands"].append(item)
            print(f"Running {script.name}; output: {log}", flush=True)
            with log.open("w") as stream:
                result = subprocess.run(
                    argv,
                    cwd=cwd,
                    env=env,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                )
            item["returncode"] = result.returncode
            if result.returncode:
                raise RuntimeError(
                    f"{script.name} failed ({result.returncode}); see {log}"
                )
        record["status"] = "commands_completed"
        print(
            f"{step.name}: commands completed. Review the generated diagnostics before proceeding."
        )
    finally:
        record["finished_utc"] = datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat()
        record_path.write_text(json.dumps(record, indent=2) + "\n")


def run_step(step: Path) -> None:
    """CLI shared by the twelve lightweight per-step entry points."""
    spec = json.loads((step / "step.json").read_text())
    parser = argparse.ArgumentParser(description=spec["title"])
    parser.add_argument(
        "--workdir",
        type=Path,
        default=Path(os.environ.get("FE2S2_WORKDIR", str(REPO / "work"))),
    )
    parser.add_argument("--action", choices=list(spec["actions"]), default="default")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only; no files or calculations",
    )
    parser.add_argument(
        "--bonds",
        nargs="+",
        type=int,
        help="Override bond dimensions for a refinement action",
    )
    args = parser.parse_args()
    work = args.workdir.expanduser().resolve()
    action = spec["actions"][args.action]
    if args.bonds and "bonds" not in action:
        parser.error("--bonds requires a refinement action")
    bonds = args.bonds or action.get("bonds", [])
    if any(bond <= 0 for bond in bonds):
        parser.error("Bond dimensions must be positive")
    commands = build_commands(action, work, bonds)
    if args.dry_run:
        print(f"Step: {step.name}; action: {args.action}; workdir: {work}")
        print("Environment:", json.dumps(spec.get("env", {}), sort_keys=True))
        print(
            "Requires:", ", ".join(action.get("requires", spec["requires"])) or "(none)"
        )
        for argv, cwd, _ in commands:
            print(f"cd {shlex.quote(str(cwd))} && {shlex.join(argv)}")
        return

    try:
        prepare_workdir(work, step, spec, action)
    except ValueError as error:
        parser.error(str(error))
    env = os.environ.copy()
    env["FE2S2_WORKDIR"] = str(work)
    env["PYTHONPATH"] = os.pathsep.join(
        path
        for path in [str(REPO / "_shared"), str(step), env.get("PYTHONPATH")]
        if path
    )
    for key, value in spec.get("env", {}).items():
        env.setdefault(key, str(value))
    execute_commands(step, args.action, commands, work, env)
