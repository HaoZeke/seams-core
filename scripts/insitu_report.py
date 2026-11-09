"""Type 7 quartiles of the live overhead samples. Writes REPORT.md."""

import json
import sys
from pathlib import Path


def type7(values, probability):
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        raise SystemExit("no samples")
    if n == 1:
        return ordered[0]
    position = (n - 1) * probability
    lo = int(position)
    frac = position - lo
    if lo + 1 >= n:
        return ordered[-1]
    return ordered[lo] * (1.0 - frac) + ordered[lo + 1] * frac


def fmt(value):
    return "%.17g" % value


def main():
    src = Path(sys.argv[1])
    dest = Path(sys.argv[2])
    data = json.loads(src.read_text())
    host = data["host"]
    lines = []
    lines.append("# In-situ overhead")
    lines.append("")
    lines.append("LAMMPS samples are its reported loop time divided by the number of steps.")
    lines.append("Every other code uses the wall clock of the same driver call, divided by the number of steps.")
    lines.append("The plain run and the in-situ run use the same input on one process and one thread.")
    lines.append("Six samples are kept. The quartiles are Hyndman-Fan type 7.")
    lines.append("Overhead is the live sample minus the plain sample at the same index, in seconds per step.")
    lines.append("")
    lines.append("Host `%s`." % host["hostname"])
    lines.append("CPU `%s`." % host["cpu"])
    lines.append("HEAD `%s`." % host["head"])
    lines.append("Tree `%s`." % host["tree"])
    lines.append("libyodaLib SHA-256 `%s`." % host["lib_sha256"])
    lines.append("LAMMPS SHA-256 `%s`." % host["lmp_sha256"])
    lines.append("GROMACS SHA-256 `%s`." % host["gmx_sha256"])
    lines.append("Plugin SHA-256 `%s`." % host.get("plugin_sha256", ""))
    lines.append("")
    lines.append("| Code | Steps | Frames | Q1 overhead | Median | Q3 | IQR |")
    lines.append("| --- | --- | --- | --- | --- | --- | --- |")
    for name, block in data["codes"].items():
        plain = block.get("plain") or []
        live = block.get("live") or []
        if len(plain) != len(live) or not plain:
            lines.append("| %s | %s | %s |  |  |  |  |" % (name, block.get("steps", ""), block.get("frames", "")))
            continue
        delta = [live[i] - plain[i] for i in range(len(plain))]
        q1 = type7(delta, 0.25)
        med = type7(delta, 0.5)
        q3 = type7(delta, 0.75)
        lines.append(
            "| %s | %s | %s | %s | %s | %s | %s |"
            % (name, block["steps"], block["frames"], fmt(q1), fmt(med), fmt(q3), fmt(q3 - q1))
        )
    lines.append("")
    for name, block in data["codes"].items():
        plain = block.get("plain") or []
        live = block.get("live") or []
        if not plain:
            continue
        lines.append("## %s" % name)
        lines.append("")
        lines.append("Plain seconds per step:")
        lines.append("")
        for value in plain:
            lines.append(fmt(value))
        lines.append("")
        lines.append("Live seconds per step:")
        lines.append("")
        for value in live:
            lines.append(fmt(value))
        lines.append("")
        lines.append("Overhead seconds per step:")
        lines.append("")
        for i in range(len(plain)):
            lines.append(fmt(live[i] - plain[i]))
        lines.append("")
    dest.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
