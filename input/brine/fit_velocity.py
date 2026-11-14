#!/usr/bin/env python3
"""Front velocity against temperature from coexistence CSVs.

Each CSV is one trajectory written by front_csv.py. The ordinary
least-squares velocity is half the slope of ice thickness against
time on every frame after --skip-ns. That standard error treats
frames as independent.

The gate uses a coarser error. The trajectory is split into
half-nanosecond windows, each window is fitted on its own, and the
reported uncertainty is the standard error of those window
velocities. The gate passes only when the cold side grows and the
warm side melts, each by more than twice that standard error.
Two seeds at the same temperature are averaged before the test.

The melting point is the zero of a weighted straight line through the
block-mean velocity at every temperature. The weight is the inverse
variance of that mean. The standard error is the delta-method error
of the zero. When the residual chi-squared per degree of freedom
exceeds one, that error is multiplied by the square root, so a line
the seeds do not support is not quoted more tightly than they allow.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path


def read_series(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("t_ns"):
            continue
        parts = line.split(",")
        rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return rows


def slope(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan"), float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    ssxx = sum((x - mx) ** 2 for x in xs)
    ssxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if ssxx <= 0:
        return float("nan"), float("nan")
    b = ssxy / ssxx
    sse = sum((y - (my + b * (x - mx))) ** 2 for x, y in zip(xs, ys))
    err = math.sqrt((sse / (n - 2)) / ssxx)
    return b, err


def temperature(path: Path) -> float:
    match = re.search(r"(\d+(?:\.\d+)?)", path.stem)
    if not match:
        raise SystemExit(f"no temperature in {path.name}")
    return float(match.group(1))


def window_velocities(rows, width):
    if not rows:
        return []
    start = rows[0][0]
    end = rows[-1][0]
    velocities = []
    left = start
    while left < end - 1e-9:
        right = min(left + width, end)
        last = abs(right - end) < 1e-9
        block = [
            row
            for row in rows
            if (left <= row[0] < right) or (last and abs(row[0] - end) < 1e-9)
        ]
        if len(block) >= 3:
            thick, _err = slope([row[0] for row in block], [row[2] for row in block])
            velocities.append(thick / 2.0)
        left = right
    return velocities


def mean_se(values):
    n = len(values)
    if n == 0:
        return float("nan"), float("nan")
    mean = sum(values) / n
    if n == 1:
        return mean, float("nan")
    var = sum((value - mean) ** 2 for value in values) / (n * (n - 1))
    return mean, math.sqrt(var)


def one_run(path: Path, skip: float, width: float):
    rows = read_series(path)
    kept = [row for row in rows if row[0] >= skip]
    thick, thick_err = slope([row[0] for row in kept], [row[2] for row in kept])
    dn, dn_err = slope([row[0] for row in kept], [row[1] for row in kept])
    block_v, block_se = mean_se(window_velocities(rows, width))
    return {
        "path": path,
        "T": temperature(path),
        "n": len(kept),
        "t_end": rows[-1][0] if rows else float("nan"),
        "v": thick / 2.0,
        "ve": thick_err / 2.0,
        "dn": dn,
        "dne": dn_err,
        "block_v": block_v,
        "block_se": block_se,
        "n_block": len(window_velocities(rows, width)),
    }


def mean_group(runs):
    n = len(runs)
    v = sum(run["v"] for run in runs) / n
    dn = sum(run["dn"] for run in runs) / n
    block_v = sum(run["block_v"] for run in runs) / n
    if n == 1:
        ve, dne = runs[0]["ve"], runs[0]["dne"]
        block_se = runs[0]["block_se"]
        kind = "slope"
    else:
        ve = math.sqrt(sum((run["v"] - v) ** 2 for run in runs) / (n * (n - 1)))
        dne = math.sqrt(sum((run["dn"] - dn) ** 2 for run in runs) / (n * (n - 1)))
        block_se = math.sqrt(
            sum((run["block_v"] - block_v) ** 2 for run in runs) / (n * (n - 1))
        )
        kind = "replicate"
    return {
        "T": runs[0]["T"],
        "nseed": n,
        "v": v,
        "ve": ve,
        "dn": dn,
        "dne": dne,
        "block_v": block_v,
        "block_se": block_se,
        "kind": kind,
        "t_end": min(run["t_end"] for run in runs),
    }


def crossing(cold, warm, key, err_key):
    dv = warm[key] - cold[key]
    if dv == 0 or math.isnan(dv):
        return float("nan"), float("nan")
    span = warm["T"] - cold["T"]
    t0 = cold["T"] - cold[key] * span / dv
    d_dv1 = -span * warm[key] / dv**2
    d_dv2 = span * cold[key] / dv**2
    ce = 0.0 if math.isnan(cold[err_key]) else cold[err_key]
    we = 0.0 if math.isnan(warm[err_key]) else warm[err_key]
    err = math.sqrt((d_dv1 * ce) ** 2 + (d_dv2 * we) ** 2)
    return t0, err


def weighted_block_zero(points):
    """Zero of the weighted line through per-temperature block means.

    Returns temperature, standard error, chi-squared, degrees of
    freedom and the factor applied to the delta-method error.
    """
    if len(points) < 2:
        return float("nan"), float("nan"), float("nan"), 0, 1.0
    sw = swt = swv = swtt = swtv = 0.0
    for point in points:
        se = point["block_se"]
        if se <= 0 or math.isnan(se):
            raise SystemExit(f"no block error at {point['T']}")
        weight = 1.0 / (se * se)
        temp = point["T"]
        vel = point["block_v"]
        sw += weight
        swt += weight * temp
        swv += weight * vel
        swtt += weight * temp * temp
        swtv += weight * temp * vel
    det = sw * swtt - swt * swt
    if det <= 0:
        return float("nan"), float("nan"), float("nan"), 0, 1.0
    intercept = (swtt * swv - swt * swtv) / det
    slope = (sw * swtv - swt * swv) / det
    if slope == 0 or math.isnan(slope):
        return float("nan"), float("nan"), float("nan"), 0, 1.0
    chi = 0.0
    for point in points:
        resid = point["block_v"] - (intercept + slope * point["T"])
        chi += (resid / point["block_se"]) ** 2
    var_a = swtt / det
    var_b = sw / det
    cov = -swt / det
    zero = -intercept / slope
    d_a = -1.0 / slope
    d_b = intercept / (slope * slope)
    var = d_a * d_a * var_a + d_b * d_b * var_b + 2.0 * d_a * d_b * cov
    err = math.sqrt(var) if var > 0 else float("nan")
    dof = len(points) - 2
    factor = 1.0
    if dof > 0 and chi > dof:
        factor = math.sqrt(chi / dof)
        err *= factor
    return zero, err, chi, dof, factor


ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "paper" / "tex" / "figures" / "data" / "gate"
MACROS = GATE / "front_macros.tex"
TEX = ROOT / "paper" / "tex" / "brine_science.tex"
# Filename stem, temperature is read from the stem, and the seed label
# is the one quoted in the table. "earlier" is the first 245 K and 255 K run.
SEED_LABELS = (
    ("pure-235.csv", "2351"),
    ("pure-235-b.csv", "2352"),
    ("pure-235-c.csv", "2353"),
    ("pure-240-a.csv", "2401"),
    ("pure-245.csv", "earlier"),
    ("pure-245-b.csv", "2452"),
    ("pure-250.csv", "2501"),
    ("pure-255.csv", "earlier"),
)
_DIGITS = (
    "Zero",
    "One",
    "Two",
    "Three",
    "Four",
    "Five",
    "Six",
    "Seven",
    "Eight",
    "Nine",
)


def gate_csvs():
    return [GATE / name for name, _label in SEED_LABELS]


def _fixed(value, digits):
    quantum = Decimal("1").scaleb(-digits)
    number = Decimal(format(value, f".{digits + 8}f"))
    return number.quantize(quantum, rounding=ROUND_HALF_UP)


def _num(value, digits):
    return f"{_fixed(value, digits):.{digits}f}"


def _pm(value, err, digits):
    shown = _fixed(value, digits)
    sign = "+" if shown >= 0 else ""
    return f"${sign}{_num(value, digits)}\\pm {_num(err, digits)}$"


def _key(temp, label):
    if label == "earlier":
        return "Earlier" + "".join(_DIGITS[int(char)] for char in str(int(temp)))
    return "".join(_DIGITS[int(char)] for char in label)


def _command(kind, key, body):
    return f"\\newcommand{{\\{kind}{key}}}{{{body}}}"


def render_macros(runs):
    by_name = {run["path"].name: run for run in runs}
    labelled = []
    for name, label in SEED_LABELS:
        run = by_name[name]
        labelled.append((run, label, _key(run["T"], label)))
    by_t = {}
    for run, _label, _key_name in labelled:
        by_t.setdefault(run["T"], []).append(run)
    points = [mean_group(by_t[temp]) for temp in sorted(by_t)]
    zero, err, chi, dof, factor = weighted_block_zero(points)
    two_235 = mean_group(
        [by_name["pure-235.csv"], by_name["pure-235-b.csv"]]
    )
    two_245 = mean_group(
        [by_name["pure-245.csv"], by_name["pure-245-b.csv"]]
    )
    one_255 = mean_group([by_name["pure-255.csv"]])
    hist_zero, hist_err = crossing(two_245, one_255, "block_v", "block_se")
    lines = [
        "% Generated by input/brine/fit_velocity.py from the gate CSVs. Do not edit.",
        _command("GateZeroTwo", "", _num(zero, 2)),
        _command("GateZeroErrTwo", "", _num(err, 2)),
        _command("GateZeroOne", "", _num(zero, 1)),
        _command("GateZeroErrOne", "", _num(err, 1)),
        _command("GateChi", "", _num(chi, 2)),
        _command("GateDof", "", str(dof)),
        _command("GateScale", "", _num(factor, 3)),
        _command("GateHistZero", "", _num(hist_zero, 1)),
        _command("GateHistErr", "", _num(hist_err, 1)),
        _command("GateMeanBlockTwo", "", _pm(two_235["block_v"], two_235["block_se"], 2)),
        _command(
            "GateMeanBlockThree",
            "",
            _pm(
                next(point["block_v"] for point in points if point["T"] == 235),
                next(point["block_se"] for point in points if point["T"] == 235),
                2,
            ),
        ),
        _command("GateMeanBlockTwoFourFive", "", _pm(two_245["block_v"], two_245["block_se"], 2)),
    ]
    for run, label, key in labelled:
        temp = str(int(run["T"])) if run["T"] == int(run["T"]) else _num(run["T"], 1)
        lines.append(_command("GateTemp", key, temp))
        lines.append(_command("GateName", key, label))
        lines.append(_command("GateOls", key, _pm(run["v"], run["ve"], 2)))
        lines.append(_command("GateBlock", key, _pm(run["block_v"], run["block_se"], 2)))
        lines.append(_command("GateDn", key, _pm(run["dn"], run["dne"], 1)))
    lines.append("")
    return "\n".join(lines) + "\n"


def macro_names(text):
    return re.findall(r"\\newcommand\{\\(\w+)\}", text)


def check_macros(text, tex):
    problems = []
    try:
        committed = MACROS.read_text()
    except FileNotFoundError:
        committed = None
        problems.append("front_macros.tex is missing")
    if committed is not None and text != committed:
        problems.append("front_macros.tex is not the bytes written from the gate CSVs")
    if "figures/data/gate/front_macros.tex" not in tex:
        problems.append("brine_science.tex does not input the gate macros")
    for name in macro_names(text):
        if "\\" + name not in tex:
            problems.append(f"brine_science.tex does not use \\{name}")
    if r"\multicolumn{3}{c}{pending}" in tex:
        problems.append("the pending 2501 row is still in the table")
    return problems


def grows(point) -> bool:
    return point["block_v"] > 2.0 * point["block_se"]


def melts(point) -> bool:
    return point["block_v"] < -2.0 * point["block_se"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csvs", nargs="*", type=Path)
    parser.add_argument("--skip-ns", type=float, default=0.06)
    parser.add_argument("--block-ns", type=float, default=0.5)
    parser.add_argument("--macros", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    paths = args.csvs or (gate_csvs() if args.macros or args.check else None)
    if not paths:
        parser.error("pass at least one CSV")
    runs = [one_run(path, args.skip_ns, args.block_ns) for path in paths]
    if args.macros or args.check:
        text = render_macros(runs)
        if args.macros:
            args.macros.parent.mkdir(parents=True, exist_ok=True)
            args.macros.write_text(text)
        if args.check:
            problems = check_macros(text, TEX.read_text())
            if problems:
                for problem in problems:
                    print(problem, file=sys.stderr)
                sys.exit(1)
            print("gate macros byte-identical")
        return
    by_t = {}
    for run in runs:
        by_t.setdefault(run["T"], []).append(run)
    points = [mean_group(by_t[t]) for t in sorted(by_t)]
    print(
        f"# skip_ns {args.skip_ns} block_ns {args.block_ns}  "
        "v is half the thickness slope, A/ns"
    )
    print(
        f"{'T':>8} {'seeds':>5} {'v':>10} {'se':>10} "
        f"{'block':>10} {'bse':>10} {'dNdt':>10} {'se':>10} {'t_end':>8}"
    )
    for point in points:
        print(
            f"{point['T']:8.1f} {point['nseed']:5d} "
            f"{point['v']:10.4f} {point['ve']:10.4f} "
            f"{point['block_v']:10.4f} {point['block_se']:10.4f} "
            f"{point['dn']:10.2f} {point['dne']:10.2f} "
            f"{point['t_end']:8.2f}"
        )
    if len(points) < 2:
        return
    cold, warm = points[0], points[-1]
    t0, t0e = crossing(cold, warm, "v", "ve")
    tb, tbe = crossing(cold, warm, "block_v", "block_se")
    cold_grows = grows(cold)
    warm_melts = melts(warm)
    print(f"# ols_zero_K {t0:.2f} {t0e:.2f}")
    print(f"# block_zero_K {tb:.2f} {tbe:.2f}")
    wz, we, chi, dof, factor = weighted_block_zero(points)
    print(
        f"# weighted_block_zero_K {wz:.2f} {we:.2f} "
        f"chi2 {chi:.2f} dof {dof} scale {factor:.3f}"
    )
    print(f"# cold_grows {str(cold_grows).lower()} warm_melts {str(warm_melts).lower()}")
    passed = cold_grows and warm_melts
    print(f"# gate {'pass' if passed else 'fail'}")
    if not passed:
        sys.exit(2)


if __name__ == "__main__":
    main()
