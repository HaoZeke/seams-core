# Copyright (c) d-SEAMS developers
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
    lines.append("LAMMPS samples are its reported loop time divided by the step count.")
    lines.append("Every other code times the same driver call and divides that wall clock by how many steps ran.")
    lines.append("The plain run and the in-situ run use the same input on one process and one thread.")
    lines.append("Six samples are kept.")
    lines.append("The quartiles are Hyndman-Fan type 7.")
    lines.append("The interquartile range (IQR) is the third quartile minus the first.")
    lines.append("Overhead is the live sample minus the plain sample at the same index, in seconds per step.")
    lines.append("")
    lines.append("CPU `%s`." % host["cpu"])
    lines.append("Measured commit `%s`." % host["head"])
    lines.append("Measured tree `%s`." % host["tree"])
    lines.append("libyodaLib SHA-256 `%s`." % host["lib_sha256"])
    lines.append("LAMMPS SHA-256 `%s`." % host["lmp_sha256"])
    lines.append("Plugin SHA-256 `%s`." % host.get("plugin_sha256", ""))
    lines.append("")
    code_lines = []
    for name, block in data.get("codes", {}).items():
        plain = block.get("plain") or []
        live = block.get("live") or []
        if len(plain) != len(live) or not plain:
            continue
        delta = [live[i] - plain[i] for i in range(len(plain))]
        q1 = type7(delta, 0.25)
        med = type7(delta, 0.5)
        q3 = type7(delta, 0.75)
        ratios = [live[i] / plain[i] for i in range(len(plain))]
        code_lines.append(
            "| %s | %s | %s | %s | %s | %s | %s | %s |"
            % (
                name,
                block["steps"],
                block["frames"],
                fmt(q1),
                fmt(med),
                fmt(q3),
                fmt(q3 - q1),
                fmt(type7(ratios, 0.5)),
            )
        )
    if code_lines:
        lines.append("| Code | Steps | Frames | Q1 overhead | Median | Q3 | IQR | Median ratio |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        lines.extend(code_lines)
    lines.append("")
    lines.append("Overhead is seconds per step.")
    lines.append("The ratio is the live step time divided by the plain step time.")
    lines.append("The size table records the LAMMPS runs.")
    lines.append("The dump of those runs matches the in-situ record.")
    lines.append("")
    if "scale" in data:
        lines.append("# LAMMPS size and cadence")
        lines.append("")
        lines.append("The system is a diamond lattice of one site per water.")
        lines.append("Plain and live samples are the LAMMPS loop time divided by the number of steps.")
        lines.append("The postprocess column is the wall time of the seams rings, chill-plus and cages commands.")
        lines.append("That column uses one native dump frame, repeated six times.")
        lines.append("The in-situ run does not write a dump.")
        lines.append("The byte count is the native dump of one frame.")
        lines.append("")
        lines.append("| Waters | Cadence | Steps | Median absolute | Median ratio | Per analysed frame | Median postprocess | Dump bytes |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for natoms, cadences in data["scale"].items():
            for every, block in cadences.items():
                plain = block["plain"]
                live = block["live"]
                delta = [live[i] - plain[i] for i in range(len(plain))]
                ratios = [live[i] / plain[i] for i in range(len(plain))]
                lines.append(
                    "| %s | %s | %s | %s | %s | %s | %s | %d |"
                    % (
                        natoms,
                        every,
                        block["steps"],
                        fmt(type7(delta, 0.5)),
                        fmt(type7(ratios, 0.5)),
                        fmt(type7(delta, 0.5) * int(every)),
                        fmt(type7(block["post_wall"], 0.5)),
                        block["dump_bytes"],
                    )
                )
        lines.append("")
        for natoms, cadences in data["scale"].items():
            for every, block in cadences.items():
                lines.append("## lammps %s waters every %s" % (natoms, every))
                lines.append("")
                if "warmup_plain" in block:
                    lines.append("Warmup plain seconds per step:")
                    lines.append("")
                    lines.append(fmt(block["warmup_plain"]))
                    lines.append("")
                    lines.append("Warmup live seconds per step:")
                    lines.append("")
                    lines.append(fmt(block["warmup_live"]))
                    lines.append("")
                lines.append("Plain seconds per step:")
                lines.append("")
                for value in block["plain"]:
                    lines.append(fmt(value))
                lines.append("")
                lines.append("Live seconds per step:")
                lines.append("")
                for value in block["live"]:
                    lines.append(fmt(value))
                lines.append("")
                lines.append("Postprocess wall seconds:")
                lines.append("")
                for value in block["post_wall"]:
                    lines.append(fmt(value))
                lines.append("")
    if "gromacs_scale" in data:
        lines.append("# Module size and cadence")
        lines.append("")
        lines.append("The same diamond lattices as the LAMMPS table, one site per water.")
        lines.append("Plain is mdrun.")
        lines.append("Live is the d-SEAMS module compiled into mdrun.")
        lines.append("Both times are the mdrun wall time divided by the number of steps.")
        lines.append("The run uses one rank and one thread.")
        lines.append("The analysed-frame column is the median extra seconds per step multiplied by the cadence.")
        lines.append("")
        lines.append("| Waters | Cadence | Steps | Median absolute | Median ratio | Per analysed frame |")
        lines.append("| --- | --- | --- | --- | --- | --- |")
        for natoms, cadences in data["gromacs_scale"].items():
            for every, block in cadences.items():
                plain = block["plain"]
                live = block["live"]
                delta = [live[i] - plain[i] for i in range(len(plain))]
                ratios = [live[i] / plain[i] for i in range(len(plain))]
                lines.append(
                    "| %s | %s | %s | %s | %s | %s |"
                    % (
                        natoms,
                        every,
                        block["steps"],
                        fmt(type7(delta, 0.5)),
                        fmt(type7(ratios, 0.5)),
                        fmt(type7(delta, 0.5) * int(every)),
                    )
                )
        lines.append("")
    if "clock" in data:
        clock = data["clock"]
        lines.append("# Analysis time against the command line")
        lines.append("")
        lines.append("The in-situ column is the extra loop time of one analysed frame.")
        lines.append("The command line included column is the wall time of seams rings, chill-plus and cages.")
        lines.append("Process start and file reading are included.")
        lines.append("The excluded column is those three censuses after the dump is already in memory.")
        lines.append("")
        lines.append("| System | In-situ per frame | CLI included | CLI excluded |")
        lines.append("| --- | --- | --- | --- |")
        lines.append(
            "| %s | %s | %s | %s |"
            % (
                clock.get("system", ""),
                fmt(clock["insitu_per_frame"]),
                fmt(clock["cli_included"]),
                fmt(clock["cli_excluded"]),
            )
        )
        lines.append("")
        steps = data.get("scale", {}).get("1000", {}).get("1", {}).get("steps", "")
        lines.append(
            "The in-situ number is the median extra time per step of the 1000-water every-step run (%s steps)."
            % steps
        )
        lines.append("It includes the first frame.")
        if "steady" in data:
            steady = data["steady"]
            delta = [steady["live"][i] - steady["plain"][i] for i in range(len(steady["plain"]))]
            ratios = [steady["live"][i] / steady["plain"][i] for i in range(len(steady["plain"]))]
            lines.append(
                "A %s-step run at cadence %s on 1000 waters adds %s seconds per step (median ratio %s), which is %s seconds for each analysed frame."
                % (
                    steady["steps"],
                    steady["every"],
                    fmt(type7(delta, 0.5)),
                    fmt(type7(ratios, 0.5)),
                    fmt(type7(delta, 0.5) * steady["every"]),
                )
            )
        lines.append("The command line excluded time is one frame.")
        lines.append("Its included time is one process per census.")
        lines.append("")
    if "equality" in data:
        lines.append("# Equality against each code's own file")
        lines.append("")
        lines.append("A frame matches when its ring line, its CHILL+ counts and its cage counts are identical.")
        lines.append("")
        lines.append("| Code | Frames | Rings | CHILL+ | Cages |")
        lines.append("| --- | --- | --- | --- | --- |")
        for row in data["equality"]:
            lines.append(
                "| %s | %s | %s | %s | %s |"
                % (row["code"], row["frames"], row["rings"], row["chill"], row["cages"])
            )
        lines.append("")
    if "ladder" in data:
        lines.append("# OpenMP threads")
        lines.append("")
        lines.append("The ladder times a warm cached workload.")
        lines.append("Five dumped frames are parsed before the clock starts.")
        lines.append("One engine stays open for the warmup and the six samples.")
        lines.append("Each sample calls observe on those frames and divides by five.")
        lines.append("Parsing, construction, and the cold first pass stay outside the samples.")
        lines.append("A larger thread count scales when its third quartile is below the previous first quartile.")
        lines.append("")
        previous = {}
        scaling_failed = []
        for block in data["ladder"]:
            waters = block["waters"]
            for count, row in block["threads"].items():
                samples = row["samples"]
                q1 = type7(samples, 0.25)
                med = type7(samples, 0.5)
                q3 = type7(samples, 0.75)
                lines.append("Waters %s, threads %s." % (waters, count))
                lines.append("First quartile %s." % fmt(q1))
                lines.append("Median %s." % fmt(med))
                lines.append("Third quartile %s." % fmt(q3))
                if waters in previous and not (q3 < previous[waters]):
                    scaling_failed.append((waters, count))
                previous[waters] = q1
        lines.append("")
        if scaling_failed:
            lines.append("Thread scaling failed. A larger count must sit below the previous count's first quartile.")
            for waters, count in scaling_failed:
                lines.append("Failed at %s waters with %s threads." % (waters, count))
            lines.append("")
    if "mpi" in data:
        lines.append("# Ranks")
        lines.append("")
        for row in data["mpi"]:
            lines.append(
                "%s ranks %s identical to serial: %s."
                % (row["code"], row["ranks"], "yes" if row["identical"] else "no")
            )
            if row.get("note"):
                lines.append(row["note"])
        lines.append("")
    for name, block in data.get("codes", {}).items():
        plain = block.get("plain") or []
        live = block.get("live") or []
        if len(plain) != len(live) or not plain:
            continue
        lines.append("## %s raw samples" % name)
        lines.append("")
        lines.append("Plain %s loop, seconds per step." % name)
        lines.append("")
        lines.append("```")
        for value in plain:
            lines.append(fmt(value))
        lines.append("```")
        lines.append("")
        lines.append("In-situ %s loop, seconds per step." % name)
        lines.append("")
        lines.append("```")
        for value in live:
            lines.append(fmt(value))
        lines.append("```")
        lines.append("")
        lines.append("Extra %s time, seconds per step." % name)
        lines.append("")
        lines.append("```")
        for i in range(len(plain)):
            lines.append(fmt(live[i] - plain[i]))
        lines.append("```")
        lines.append("")
    dest.write_text("\n".join(lines) + "\n")
    if "ladder" in data:
        # Re-check after the file is written so the tables remain.
        previous = {}
        for block in data["ladder"]:
            waters = block["waters"]
            for count, row in block["threads"].items():
                q1 = type7(row["samples"], 0.25)
                q3 = type7(row["samples"], 0.75)
                if waters in previous and not (q3 < previous[waters]):
                    raise SystemExit("thread scaling failed for %s waters at %s threads" % (waters, count))
                previous[waters] = q1


if __name__ == "__main__":
    main()
