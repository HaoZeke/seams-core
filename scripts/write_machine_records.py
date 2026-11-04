#!/usr/bin/env python3
"""Write the machine records the paper statistics read.

The processor name is stored without a register mark, as
'Intel Xeon Processor'. Run this before paper/tex/figures/stats.py
when those records are not already on disk.
"""

import json
import os
import platform
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "paper" / "tex" / "figures" / "data"


def published_cpu_model(name: str) -> str:
    return " ".join(name.replace("(R)", "").split())


def cpuinfo() -> dict[str, str]:
    fields: dict[str, str] = {}
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        fields.setdefault(key.strip(), value.strip())
    return fields


def host_text() -> str:
    info = cpuinfo()
    gcc = subprocess.run(
        ["c++", "--version"], check=False, capture_output=True, text=True
    ).stdout.splitlines()
    gcc_line = gcc[0] if gcc else "c++ (GCC) unknown"
    lines = [
        f"uname -m: {platform.machine()}",
        f"nproc: {os.cpu_count() or 1}",
        f"model name\t: {published_cpu_model(info.get('model name', 'unknown'))}",
        f"cpu family\t: {info.get('cpu family', '')}",
        f"model\t\t: {info.get('model', '')}",
        gcc_line,
    ]
    return "\n".join(lines) + "\n"


def environment() -> dict:
    lscpu = subprocess.run(
        ["lscpu"], check=False, capture_output=True, text=True
    ).stdout.replace("(R)", "")
    uname = {
        key: value.replace("(R)", "") if isinstance(value, str) else value
        for key, value in platform.uname()._asdict().items()
        if key != "node"
    }
    return {
        "uname": uname,
        "lscpu": lscpu,
    }


def main() -> None:
    text = host_text()
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "host.txt").write_text(text)
    (DATA / "churn-host.txt").write_text(text)
    (DATA / "h2h").mkdir(parents=True, exist_ok=True)
    (DATA / "h2h" / "environment.json").write_text(
        json.dumps(environment(), indent=2) + "\n"
    )
    print("wrote host.txt, churn-host.txt and h2h/environment.json")


if __name__ == "__main__":
    main()
