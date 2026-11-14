#!/usr/bin/env python3
"""Front velocity against temperature for the pure-water coexistence CSVs.

The ordinary-least-squares velocity and the half-nanosecond blocks are
the same quantities as input/brine/fit_velocity.py: skip the first
0.06 ns, take half the thickness slope, and split the whole trajectory
into 0.5 ns windows. The OLS standard error treats frames as
independent.

The autocorrelation-corrected standard error is a Bartlett-kernel
Newey-West estimator of that same slope. The lag is the Sokal window
of the residual autocorrelation, the smallest k with k greater than
six times the running integrated autocorrelation time. The score of
the demeaned regression is x*e, and the kernel weight at lag k is
1 - k/(L+1). The reported velocity error is half of that slope error,
because the front velocity is half the thickness slope. dN/dt uses the
same estimator on the crystal count.

One point is drawn per temperature. Two seeds are averaged. The error
bar is the autocorrelation-corrected standard error of that mean.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_fit():
    path = ROOT / "input" / "brine" / "fit_velocity.py"
    spec = importlib.util.spec_from_file_location("fit_velocity", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def integrated_time(resid):
    """Sokal window. Returns (tau_int, lag)."""
    n = len(resid)
    g0 = sum(e * e for e in resid)
    if n < 4 or g0 <= 0:
        return 0.5, 1

    def rho(k):
        return sum(resid[i] * resid[i - k] for i in range(k, n)) / g0

    acc = 0.5
    for k in range(1, n // 2):
        acc += rho(k)
        if k >= 6.0 * max(acc, 0.5):
            return max(acc, 0.5), k
    return max(acc, 0.5), max(1, n // 2 - 1)


def hac_slope(xs, ys):
    """Return (slope, naive SE, Newey-West SE, lag, tau)."""
    n = len(xs)
    if n < 4:
        return float("nan"), float("nan"), float("nan"), 0, float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    xc = [x - mx for x in xs]
    ssxx = sum(x * x for x in xc)
    if ssxx <= 0:
        return float("nan"), float("nan"), float("nan"), 0, float("nan")
    slope = sum(x * (y - my) for x, y in zip(xc, ys)) / ssxx
    resid = [(y - my) - slope * x for x, y in zip(xc, ys)]
    sse = sum(e * e for e in resid)
    naive = math.sqrt((sse / (n - 2)) / ssxx)
    tau, lag = integrated_time(resid)
    score = [x * e for x, e in zip(xc, resid)]
    kernel = sum(v * v for v in score)
    for k in range(1, lag + 1):
        weight = 1.0 - k / (lag + 1.0)
        gamma = sum(score[i] * score[i - k] for i in range(k, n))
        kernel += 2.0 * weight * gamma
    hac = math.sqrt(max(kernel, 0.0) / ssxx**2)
    return slope, naive, hac, lag, tau


def analyse(fit, path: Path, skip: float, width: float):
    rows = fit.read_series(path)
    kept = [row for row in rows if row[0] >= skip]
    thick, _naive, hac, lag, tau = hac_slope(
        [row[0] for row in kept], [row[2] for row in kept]
    )
    dn, _dn_naive, dn_hac, _dn_lag, _dn_tau = hac_slope(
        [row[0] for row in kept], [row[1] for row in kept]
    )
    base = fit.one_run(path, skip, width)
    base["hac"] = hac / 2.0
    base["lag"] = lag
    base["tau"] = tau
    base["dn_hac"] = dn_hac
    # one_run already stored v from the same slope. Keep that value.
    return base


def group(runs):
    n = len(runs)
    v = sum(run["v"] for run in runs) / n
    dn = sum(run["dn"] for run in runs) / n
    block = sum(run["block_v"] for run in runs) / n
    # Independent seeds: the mean's variance is the mean of the variances.
    hac = math.sqrt(sum(run["hac"] ** 2 for run in runs) / n**2)
    dn_hac = math.sqrt(sum(run["dn_hac"] ** 2 for run in runs) / n**2)
    if n == 1:
        block_se = runs[0]["block_se"]
    else:
        block_se = math.sqrt(
            sum((run["block_v"] - block) ** 2 for run in runs) / (n * (n - 1))
        )
    return {
        "T": runs[0]["T"],
        "nseed": n,
        "v": v,
        "hac": hac,
        "dn": dn,
        "dn_hac": dn_hac,
        "block": block,
        "block_se": block_se,
        "n_block": min(run["n_block"] for run in runs),
        "t_end": min(run["t_end"] for run in runs),
    }


def pdf_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def write_pdf(points, path: Path):
    """One panel, Helvetica, error bars on the OLS velocity."""
    width, height = 460, 320
    left, right, bottom, top = 58, 18, 46, 28
    plot_w = width - left - right
    plot_h = height - bottom - top
    temps = [point["T"] for point in points]
    t0, t1 = min(temps) - 4.0, max(temps) + 4.0
    velocities = [point["v"] for point in points]
    errors = [point["hac"] for point in points]
    low = min(v - e for v, e in zip(velocities, errors))
    high = max(v + e for v, e in zip(velocities, errors))
    pad = 0.15 * max(high - low, 0.5)
    y0, y1 = min(low, 0.0) - pad, max(high, 0.0) + pad

    def xp(t):
        return left + (t - t0) / (t1 - t0) * plot_w

    def yp(v):
        return bottom + (v - y0) / (y1 - y0) * plot_h

    commands = []
    # axes
    commands.append("0.75 w")
    commands.append(f"{left} {bottom} m {left} {height - top} l S")
    commands.append(f"{left} {bottom} m {width - right} {bottom} l S")
    zero = yp(0.0)
    commands.append("0.6 w [2 2] 0 d")
    commands.append(f"{left} {zero:.2f} m {width - right} {zero:.2f} l S")
    commands.append("[] 0 d")
    # ticks
    commands.append("0.6 w")
    t = math.ceil(t0 / 5.0) * 5
    while t <= t1 + 1e-9:
        x = xp(t)
        commands.append(f"{x:.2f} {bottom} m {x:.2f} {bottom - 4} l S")
        commands.append(
            f"BT /F1 9 Tf {x - 8:.2f} {bottom - 16:.2f} Td ({int(t)}) Tj ET"
        )
        t += 5
    y_step = 1.0 if (y1 - y0) > 4 else 0.5
    y = math.ceil(y0 / y_step) * y_step
    while y <= y1 + 1e-9:
        yy = yp(y)
        commands.append(f"{left} {yy:.2f} m {left - 4} {yy:.2f} l S")
        label = f"{y:.1f}"
        commands.append(
            f"BT /F1 9 Tf {left - 32:.2f} {yy - 3:.2f} Td ({pdf_escape(label)}) Tj ET"
        )
        y += y_step
    commands.append(
        f"BT /F1 11 Tf {left + plot_w / 2 - 48:.2f} 16 Td (Temperature (K)) Tj ET"
    )
    # 90 degrees counter-clockwise, baseline starting at (16, 78).
    commands.append(
        "BT /F1 11 Tf 0 1 -1 0 16 78 Tm (Front velocity (A/ns)) Tj ET"
    )
    # points and error bars, Okabe-Ito blue
    commands.append("0 0.447 0.698 RG 0 0.447 0.698 rg")
    commands.append("1.1 w")
    for point in points:
        x = xp(point["T"])
        y = yp(point["v"])
        ylo = yp(point["v"] - point["hac"])
        yhi = yp(point["v"] + point["hac"])
        commands.append(f"{x:.2f} {ylo:.2f} m {x:.2f} {yhi:.2f} l S")
        commands.append(f"{x - 4:.2f} {ylo:.2f} m {x + 4:.2f} {ylo:.2f} l S")
        commands.append(f"{x - 4:.2f} {yhi:.2f} m {x + 4:.2f} {yhi:.2f} l S")
        commands.append(f"{x - 2.4:.2f} {y - 2.4:.2f} 4.8 4.8 re f")
    stream = "\n".join(commands).encode()
    objects = []

    def add(body: bytes):
        objects.append(body)

    add(b"<< /Type /Catalog /Pages 2 0 R >>")
    add(
        f"<< /Type /Pages /Count 1 /Kids [3 0 R] >>".encode()
    )
    add(
        (
            f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {width} {height}] "
            f"/Contents 4 0 R /Resources << /Font << /F1 5 0 R >> >> >>"
        ).encode()
    )
    add(b"<< /Length %d >>\nstream\n" % len(stream) + stream + b"\nendstream")
    add(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    out = bytearray(b"%PDF-1.4\n")
    offsets = [0]
    for i, body in enumerate(objects, start=1):
        offsets.append(len(out))
        out += f"{i} 0 obj\n".encode() + body + b"\nendobj\n"
    xref = len(out)
    out += f"xref\n0 {len(objects) + 1}\n".encode()
    out += b"0000000000 65535 f \n"
    for offset in offsets[1:]:
        out += f"{offset:010d} 00000 n \n".encode()
    out += (
        f"trailer << /Size {len(objects) + 1} /Root 1 0 R >>\n"
        f"startxref\n{xref}\n%%EOF\n"
    ).encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csvs", nargs="+", type=Path)
    parser.add_argument("--skip-ns", type=float, default=0.06)
    parser.add_argument("--block-ns", type=float, default=0.5)
    parser.add_argument("--figure", type=Path, default=None)
    args = parser.parse_args()
    fit = load_fit()
    runs = [analyse(fit, path, args.skip_ns, args.block_ns) for path in args.csvs]
    by_t = {}
    for run in runs:
        by_t.setdefault(run["T"], []).append(run)
    points = [group(by_t[t]) for t in sorted(by_t)]
    print(
        f"# skip_ns {args.skip_ns} block_ns {args.block_ns}  "
        "v is half the thickness slope, A/ns"
    )
    print(
        f"# hac is the Newey-West standard error of v; "
        "lag is the Sokal window of the thickness residuals"
    )
    print(
        f"{'file':<16} {'T':>6} {'v':>8} {'ols_se':>8} {'hac_se':>8} "
        f"{'block':>8} {'bse':>8} {'n':>4} {'dNdt':>8} {'dN_hac':>8} {'t':>6}"
    )
    for run in runs:
        print(
            f"{run['path'].name:<16} {run['T']:6.1f} "
            f"{run['v']:8.4f} {run['ve']:8.4f} {run['hac']:8.4f} "
            f"{run['block_v']:8.4f} {run['block_se']:8.4f} {run['n_block']:4d} "
            f"{run['dn']:8.2f} {run['dn_hac']:8.2f} {run['t_end']:6.2f} "
            f"lag {run['lag']} tau {run['tau']:.2f}"
        )
    print("# one point per temperature; hac_se is the standard error of the mean")
    print(
        f"{'T':>6} {'seeds':>5} {'v':>8} {'hac_se':>8} {'block':>8} {'bse':>8} "
        f"{'dNdt':>8} {'dN_hac':>8} {'t':>6}"
    )
    for point in points:
        print(
            f"{point['T']:6.1f} {point['nseed']:5d} "
            f"{point['v']:8.4f} {point['hac']:8.4f} "
            f"{point['block']:8.4f} {point['block_se']:8.4f} "
            f"{point['dn']:8.2f} {point['dn_hac']:8.2f} {point['t_end']:6.2f}"
        )
    if args.figure:
        write_pdf(points, args.figure)
        print(args.figure, args.figure.stat().st_size)


if __name__ == "__main__":
    main()
