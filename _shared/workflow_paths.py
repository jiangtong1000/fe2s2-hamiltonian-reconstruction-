from pathlib import Path
import os

CODE_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = (
    Path(os.environ.get("FE2S2_WORKDIR", str(CODE_ROOT / "work")))
    .expanduser()
    .resolve()
)
