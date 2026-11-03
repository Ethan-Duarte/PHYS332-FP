"""Launch the local UI, keeping any Python startup traceback in a text file."""
from pathlib import Path
import runpy
import traceback

root = Path(__file__).resolve().parent
try:
    runpy.run_path(str(root / "app.py"), run_name="__main__")
except BaseException:
    error = traceback.format_exc()
    (root / "startup_error.log").write_text(error, encoding="utf-8")
    print("\nN-Baddies could not start. Error saved to startup_error.log:\n")
    print(error, flush=True)
    raise SystemExit(1)
