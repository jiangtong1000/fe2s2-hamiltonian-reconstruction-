"""Meaningful packaging checks; no SCF/DMRG/factorization imports or execution.

Tests that initialize a work directory need the external reference FCIDUMP. They use
FE2S2_REFERENCE_FCIDUMP when it is set; otherwise the file is downloaded once from its
pinned upstream commit, and those tests are skipped if that is not possible.
"""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from check_repository import verify_structure

SPEC = json.loads((ROOT / "01_inputs/inputs/external_inputs.json").read_text())[
    "reference.FCIDUMP"
]
CACHE = None
FETCH_ERROR = None


def setUpModule():
    """Provide one verified local copy of the reference FCIDUMP to every step-01 run."""
    global CACHE, FETCH_ERROR
    if os.environ.get(SPEC["local_copy_variable"]):
        return
    sys.path.insert(0, str(ROOT / "01_inputs"))
    from prepare_inputs import reference_fcidump

    CACHE = Path(tempfile.mkdtemp(prefix="fe2s2 reference "))
    try:
        reference_fcidump(CACHE / "fe2s2")
    except OSError as error:
        FETCH_ERROR = f"reference FCIDUMP unavailable ({error})"
        return
    os.environ[SPEC["local_copy_variable"]] = str(CACHE / "fe2s2")


def tearDownModule():
    if CACHE is not None:
        os.environ.pop(SPEC["local_copy_variable"], None)
        shutil.rmtree(CACHE, ignore_errors=True)


def needs_reference(test):
    def wrapper(self):
        if FETCH_ERROR:
            self.skipTest(FETCH_ERROR)
        return test(self)

    wrapper.__name__ = test.__name__
    return wrapper


class PackagingTests(unittest.TestCase):
    def call(self, step, *args):
        return subprocess.run(
            [sys.executable, str(ROOT / step / "run.py"), *map(str, args)],
            cwd="/tmp",
            text=True,
            capture_output=True,
        )

    def test_structure_and_copied_hashes(self):
        self.assertEqual(verify_structure(), 12)

    def test_every_action_dry_run_is_side_effect_free(self):
        with tempfile.TemporaryDirectory(prefix="fe2s2 packaging ") as folder:
            work = Path(folder) / "never created"
            for step in sorted(ROOT.glob("[0-9][0-9]_*")):
                spec = json.loads((step / "step.json").read_text())
                for name, action in spec["actions"].items():
                    with self.subTest(step=step.name, action=name):
                        args = ["--workdir", work, "--action", name, "--dry-run"]
                        if "bonds" in action:
                            args += ["--bonds", "123", "234"]
                        result = self.call(step.name, *args)
                        self.assertEqual(result.returncode, 0, result.stderr)
                        self.assertIn(str(work), result.stdout)
                        for command in action["commands"]:
                            self.assertIn(command["script"], result.stdout)
                        if "bonds" in action:
                            self.assertIn("--refine-bonds 123 234", result.stdout)
                        self.assertFalse(work.exists())

    @needs_reference
    def test_initialization_and_overwrite_protection(self):
        with tempfile.TemporaryDirectory(prefix="fe2s2 packaging ") as folder:
            work = Path(folder) / "new run"
            result = self.call("01_inputs", "--workdir", work)
            self.assertEqual(result.returncode, 0, result.stderr)
            ref = work / "reference_zhendong/FCIDUMP"
            before = hashlib.sha256(ref.read_bytes()).hexdigest()
            self.assertEqual(before, SPEC["sha256"])
            provenance = json.loads(
                (work / "reference_zhendong/provenance.json").read_text()
            )
            self.assertEqual(before, provenance["fcidump_sha256"])
            self.assertEqual(provenance["scalar_in_fcidump"], 0.0)
            self.assertEqual(provenance["external_core_energy"], -4976.26532397)
            self.assertTrue((work / "inputs/basis.json").exists())
            records = list((work / "run_logs").glob("*.json"))
            self.assertEqual(
                json.loads(records[0].read_text())["status"], "commands_completed"
            )
            rerun = self.call("01_inputs", "--workdir", work)
            self.assertNotEqual(rerun.returncode, 0)
            self.assertIn("Outputs already exist", rerun.stderr)
            self.assertEqual(hashlib.sha256(ref.read_bytes()).hexdigest(), before)

    def test_reference_with_wrong_checksum_is_refused(self):
        with tempfile.TemporaryDirectory(prefix="fe2s2 packaging ") as folder:
            fake = Path(folder) / "fe2s2"
            fake.write_text("not the reference FCIDUMP\n")
            work = Path(folder) / "run"
            result = subprocess.run(
                [sys.executable, str(ROOT / "01_inputs/run.py"), "--workdir", work],
                cwd="/tmp",
                text=True,
                capture_output=True,
                env=dict(os.environ, **{SPEC["local_copy_variable"]: str(fake)}),
            )
            self.assertNotEqual(result.returncode, 0)
            record = json.loads(next((work / "run_logs").glob("*.json")).read_text())
            self.assertEqual(record["status"], "failed")
            log = (work / record["commands"][0]["log"]).read_text()
            self.assertIn("SHA-256 mismatch", log)
            self.assertFalse((work / "reference_zhendong/FCIDUMP").exists())
            self.assertFalse((work / "geometry.xyz").exists())

    def test_unmarked_nonempty_directory_is_untouched(self):
        with tempfile.TemporaryDirectory() as folder:
            work = Path(folder)
            sentinel = work / "original_result.txt"
            sentinel.write_text("keep")
            result = self.call("01_inputs", "--workdir", work)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("nonempty unmarked", result.stderr)
            self.assertEqual(list(work.iterdir()), [sentinel])
            self.assertEqual(sentinel.read_text(), "keep")

    @needs_reference
    def test_missing_prerequisite_fails_before_scientific_imports(self):
        with tempfile.TemporaryDirectory() as folder:
            work = Path(folder) / "run"
            self.assertEqual(self.call("01_inputs", "--workdir", work).returncode, 0)
            result = self.call("03_localize_select", "--workdir", work)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Missing prerequisite", result.stderr)
            self.assertNotIn("ModuleNotFoundError", result.stderr)

    def test_invalid_refinement_options_do_not_create_files(self):
        with tempfile.TemporaryDirectory() as folder:
            work = Path(folder) / "unused"
            for args, message in [
                (["--bonds", "2000"], "requires a refinement action"),
                (["--action", "refine-cas20", "--bonds", "0"], "must be positive"),
            ]:
                result = self.call("09_production_states", "--workdir", work, *args)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stderr)
                self.assertFalse(work.exists())

    def test_failed_command_is_logged_and_stops_the_action(self):
        sys.path.insert(0, str(ROOT / "_shared"))
        from step_runner import execute_commands

        with tempfile.TemporaryDirectory() as folder:
            work = Path(folder)
            (work / "run_logs").mkdir()
            script = work / "fail.py"
            script.write_text(
                "print('diagnostic before failure'); raise SystemExit(7)\n"
            )
            commands = [
                ([sys.executable, str(script)], work, script),
                ([sys.executable, "-c", "raise SystemExit(99)"], work, script),
            ]
            with self.assertRaisesRegex(RuntimeError, r"failed \(7\)"):
                execute_commands(ROOT / "01_inputs", "test-failure", commands, work, {})
            record = json.loads(next((work / "run_logs").glob("*.json")).read_text())
            self.assertEqual(record["status"], "failed")
            self.assertEqual(len(record["commands"]), 1)
            self.assertEqual(record["commands"][0]["returncode"], 7)
            self.assertIn(
                "diagnostic before failure",
                (work / record["commands"][0]["log"]).read_text(),
            )
            self.assertIn("finished_utc", record)

    def test_semantic_partition_retained(self):
        maps = json.loads(
            (ROOT / "06_natural_orbitals/results/orbitals/index_maps.json").read_text()
        )
        self.assertEqual(
            [len(maps[k]) for k in ["I", "K6", "A14", "V"]], [34, 6, 14, 564]
        )
        self.assertEqual(
            sorted(sum([maps[k] for k in ["I", "K6", "A14", "V"]], [])),
            list(range(618)),
        )
        self.assertEqual(maps["nelec_full"], 98)
        self.assertFalse(maps["fixed_trial_occupations_restrict_full_propagation"])
        embeddings = json.loads(
            (
                ROOT / "06_natural_orbitals/results/orbitals/trial_embeddings.json"
            ).read_text()
        )
        self.assertEqual(len(embeddings["cas20"]["fixed_double_retained"]), 34)
        self.assertEqual(len(embeddings["cas14"]["fixed_double_retained"]), 40)
        self.assertTrue(embeddings["cas14"]["occupations_are_trial_properties_only"])


if __name__ == "__main__":
    unittest.main()
