"""Hardware, compiler and library versions for the serial campaign."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import ENV, RECORDS, ROOT, loadavg, scrub

COMMANDS = [
    ["uname", "-srm"],
    ["lscpu"],
    ["awk", "/MemTotal|MemFree|MemAvailable|SwapTotal/ {print}", "/proc/meminfo"],
    ["g++", "--version"],
    ["meson", "--version"],
    ["pkg-config", "--modversion", "eigen3"],
    ["pkg-config", "--modversion", "flexiblas"],
    ["git", "rev-parse", "HEAD"],
    ["git", "-C", "../seams-chill-published", "rev-parse", "HEAD"],
]


def main() -> None:
    RECORDS.mkdir(parents=True, exist_ok=True)
    chunks = [f"loadavg {loadavg()}"]
    for cmd in COMMANDS:
        proc = subprocess.run(cmd, cwd=ROOT, env=ENV, text=True, capture_output=True)
        chunks.append("command " + " ".join(cmd))
        chunks.append(f"returncode {proc.returncode}")
        chunks.append(proc.stdout)
        chunks.append(proc.stderr)
    py = ROOT / "benchmarks" / "cpc" / ".venv" / "bin" / "python"
    code = (
        "import numpy,ovito\n"
        "print('numpy', numpy.__version__)\n"
        "print('ovito', '.'.join(map(str, ovito.version)))\n"
    )
    proc = subprocess.run([str(py), "-c", code], cwd=ROOT, env=ENV, text=True, capture_output=True)
    chunks.append("command python -c versions")
    chunks.append(proc.stdout)
    chunks.append(proc.stderr)
    text = scrub("\n".join(chunks))
    (RECORDS / "environment.txt").write_text(text)
    print("environment.txt")


if __name__ == "__main__":
    main()
