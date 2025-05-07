from morphoparse import cli, core
from pathlib import Path
import sys

def main():
    try:
        args = cli.get_args()
        core.run_pipeline(args)
    except Exception as e:
        try:
            log_dir = Path(getattr(args, "output", "") or getattr(args, "input", "")).resolve().parent
        except Exception:
            log_dir = Path.cwd()
        print(f"[ERROR] {e}", file=sys.stderr)
        with open(log_dir / "morphoparse_errors.log", "a") as log:
            log.write(f"{e}\n")
    
    sys.exit(0)

if __name__ == "__main__":
    main()