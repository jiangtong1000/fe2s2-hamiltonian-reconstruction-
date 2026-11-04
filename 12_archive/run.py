from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_shared"))
from step_runner import run_step

if __name__ == "__main__":
    run_step(Path(__file__).resolve().parent)
