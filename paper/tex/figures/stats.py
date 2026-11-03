"""Rounded LaTeX macros and table rows from the measured records.

Reads paper/tex/figures/data/ and writes data/macros.tex, the table bodies
data/tab-*.tex and the plot tables data/plot-*.dat. The fallback thresholds
come from their defaults in src/. Standard library only.

Medians carry exact binomial order-statistic intervals at 95% coverage and
quartiles interpolated linearly between order statistics. Ratios of medians
carry percentile intervals from 10,000 bootstrap draws with seed 20261103,
the convention of the head-to-head campaign.

Usage: python3 stats.py [DATA_DIR]
"""

import csv
import decimal
import json
import math
import pathlib
import random
import re
import statistics
import sys

DATA = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else
                    pathlib.Path(__file__).parent / "data")
H2H = DATA / "h2h"
ROOT = pathlib.Path(__file__).resolve().parents[3]
SOURCE = ROOT / "src"
SEED = 20261103
DRAWS = 10000
COVERAGE = 0.95
TAIL = (1.0 - COVERAGE) / 2.0
MS_PER_S = 1000.0
# The headline rung of the version comparison, whose numbers the abstract,
# the summary and the conclusions print with HEADLINE_FIGURES figures.
HEADLINE_RUNG = 1000
HEADLINE_FIGURES = 2
# Ideal ice Ic has two six-rings per molecule: each molecule lies on twelve
# and each ring has six.
SIX_RINGS_PER_MOLECULE = 2
# The ring columns of the records count the rings of MIN_RING to MAX_RING
# members, and the prism columns the blocks of the square, pentagonal and
# hexagonal tubes.
MIN_RING = 3
MAX_RING = 6
RING_SIZES = range(MIN_RING, MAX_RING + 1)
PRISM_SIZES = (4, 5, 6)
# cells() fills CELL_COLUMNS table columns, a sample with its call count
# SAMPLE_COLUMNS and a ratio of medians with its interval RATIO_COLUMNS.
CELL_COLUMNS = 3
SAMPLE_COLUMNS = 1 + CELL_COLUMNS
RATIO_COLUMNS = 2
# The sweep record prints the medians with three decimals and the calls with
# four, so a median of the calls lies within half a unit of the third decimal
# plus half a unit of the fourth of the printed one.
ROW_TOLERANCE_MS = 0.00055
# A p-value below PVALUE_FIXED prints in scientific notation.
PVALUE_FIXED = 1e-3
# A mark of 1.8 pt spans a factor of MARK_SPAN["churn"] on the axis of the
# update figure and more than MARK_SPAN["other"] on the axes of the others.
MARK_SPAN = {"churn": 1.1, "other": 1.3}
# The text gives counts below eleven in words.
WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight",
         "nine", "ten")
# The captions of these figures give the fewest calls behind one of their marks.
CALL_FIGURES = ("versions", "codes")

macros = []
# The calls behind each mark of a plot file, in the order of its rows.
plot_calls = {}


class Rounded(str):
    """A rounded number that remembers the value it was rounded from."""

    def __new__(cls, text, value):
        obj = super().__new__(cls, text)
        obj.value = value
        return obj


def sig(x, figures=3):
    if x == 0:
        return Rounded("0", x)
    digits = figures - 1 - int(math.floor(math.log10(abs(x))))
    value = round(x, digits)
    if digits <= 0:
        return Rounded(str(int(round(value))), x)
    return Rounded(f"{value:.{digits}f}", x)


def half_up(text, figures):
    """The rounding of siunitx (round-half = up) applied to a printed number."""
    d = decimal.Decimal(text)
    if d == 0:
        return d
    step = decimal.Decimal(1).scaleb(d.adjusted() - figures + 1)
    return d.quantize(step, rounding=decimal.ROUND_HALF_UP)


def define(name, body):
    macros.append(f"\\expandafter\\def\\csname dsv@{name}\\endcsname{{{body}}}")
    # \rtwo and \rthree round the printed string again; where that differs
    # from rounding the measured value, the preamble reads these instead.
    if isinstance(body, Rounded):
        for figures, suffix in ((2, "sig2"), (3, "sig3")):
            direct = sig(body.value, figures)
            if half_up(body, figures) != decimal.Decimal(direct):
                macros.append(f"\\expandafter\\def\\csname dsv@{name}@{suffix}"
                              f"\\endcsname{{{direct}}}")


def word(n):
    return WORDS[int(n)]


def median(xs):
    s = sorted(xs)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else 0.5 * (s[mid - 1] + s[mid])


def median_interval(xs):
    """Largest k with P(Bin(n,1/2) < k) <= (1-coverage)/2; [x_(k), x_(n+1-k)]."""
    s = sorted(xs)
    n = len(s)
    tail = 0.0
    k = 0
    while k < n:
        nxt = tail + math.comb(n, k) / 2.0 ** n
        if nxt > TAIL:
            break
        tail = nxt
        k += 1
    k = max(k, 1)
    return s[k - 1], s[n - k]


def quartiles(xs):
    """First and third quartiles, interpolated linearly between order statistics."""
    q1, _, q3 = statistics.quantiles(xs, n=4, method="inclusive")
    return q1, q3


def span(q1, q3):
    return f"\\iqr{{{sig(q1)}}}{{{sig(q3)}}}"


def cells(xs):
    """Median, 95% interval and quartiles of one sample, three table cells."""
    lo, hi = median_interval(xs)
    return (f"\\num{{{sig(median(xs))}}} & [\\num{{{sig(lo)}}}, "
            f"\\num{{{sig(hi)}}}] & {span(*quartiles(xs))}")


def not_run(columns, mark="a"):
    return f"\\multicolumn{{{columns}}}{{c}}{{not run\\textsuperscript{{{mark}}}}}"


def write_table(name, rows):
    """A table body, between the markers that keep snapper from joining its rows."""
    (DATA / name).write_text("\n".join(["% snapper:off", *rows, "% snapper:on"]) + "\n")


def plot_row(x, xs):
    lo, hi = median_interval(xs)
    return f"{x} {median(xs)} {lo} {hi} " + " ".join(str(q) for q in quartiles(xs))


def ratio_interval(num, den, rng):
    ratios = []
    for _ in range(DRAWS):
        a = [num[rng.randrange(len(num))] for _ in num]
        b = [den[rng.randrange(len(den))] for _ in den]
        ratios.append(median(a) / median(b))
    ratios.sort()
    lo = ratios[int(math.floor(TAIL * DRAWS))]
    hi = ratios[int(math.ceil((1.0 - TAIL) * DRAWS)) - 1]
    return lo, hi


def read_samples(paths):
    columns = None
    rows = []
    for path in paths:
        for line in path.read_text().splitlines():
            if line.startswith("columns "):
                columns = line.split()[1:]
                continue
            parts = line.split()
            if parts and all(_is_float(p) for p in parts) and columns and \
                    len(parts) == len(columns):
                rows.append([float(p) for p in parts])
    return columns, rows


def _is_float(text):
    try:
        float(text)
    except ValueError:
        return False
    return True


def relative(m, lo, hi):
    """Half-width of the larger side of an interval, percent of the centre."""
    return sig(100.0 * max(hi - m, m - lo) / m, 2)


def define_median(base, xs):
    m = median(xs)
    lo, hi = median_interval(xs)
    q1, q3 = quartiles(xs)
    define(base, sig(m))
    define(base + "@lo", sig(lo))
    define(base + "@hi", sig(hi))
    define(base + "@rel", relative(m, lo, hi))
    define(base + "@q1", sig(q1))
    define(base + "@q3", sig(q3))
    define(base + "@n", str(len(xs)))
    return m, lo, hi


def slope(points):
    """Least-squares slope of log(time) on log(sites)."""
    xs = [math.log(n) for n, _ in points]
    ys = [math.log(t) for _, t in points]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def read_counts(path):
    parts = path.read_text().split()
    return dict(zip(parts[0::2], parts[1::2]))


STAGES = ("ring", "hc", "ddc", "path")
PATH = STAGES.index("path")


def versions():
    """d-SEAMS 1.0.1 against 2.x on ideal ice Ic frames, per-call samples."""
    rng = random.Random(SEED)
    sizes2 = sorted(int(p.stem.split("-")[-1])
                    for p in DATA.glob("dseams2-ic-*.txt"))
    sizes1 = sorted(int(p.stem.split("-")[-1])
                    for p in DATA.glob("dseams1-ic-*.txt"))
    samples = {}
    for tag, sizes in (("one", sizes1), ("two", sizes2)):
        version = "dseams1" if tag == "one" else "dseams2"
        for n in sizes:
            cols, rows = read_samples([DATA / f"{version}-ic-{n}.txt"])
            assert cols == ["ring_ms", "hc_ms", "ddc_ms", "path_ms"], cols
            samples[(tag, n)] = rows
            counts = read_counts(DATA / f"{version}-ic-{n}.err")
            for key, value in counts.items():
                define(f"v@{n}@count@{tag}@{key}", f"\\num{{{value}}}")
            for j, stage in enumerate(STAGES):
                define_median(f"v@{n}@{stage}@{tag}", [r[j] for r in rows])
    common = [n for n in sizes1 if n in sizes2]
    define("v@largest", f"\\num{{{common[-1]}}}")
    define("v@smallest", f"\\num{{{common[0]}}}")
    define("v@largest@two", f"\\num{{{sizes2[-1]}}}")
    define("v@rungs@one", str(len(sizes1)))
    define("v@rungs@two", str(len(sizes2)))
    minimum = min(len(samples[("one", n)]) for n in sizes1)
    define("v@minreps@one", str(minimum))
    lines = []
    bounds = {}
    for n in common:
        one = samples[("one", n)]
        two = samples[("two", n)]
        for j, stage in enumerate(STAGES):
            a = [r[j] for r in one]
            b = [r[j] for r in two]
            lo, hi = ratio_interval(a, b, rng)
            bounds[(n, stage)] = (lo, hi)
            r = median(a) / median(b)
            base = f"v@{n}@{stage}@ratio"
            define(base, sig(r))
            define(base + "@lo", sig(lo))
            define(base + "@hi", sig(hi))
            define(base + "@rel", relative(r, lo, hi))
            define(f"v@{n}@{stage}@separated",
                   "yes" if min(a) > max(b) else "no")
    # The caption of the version table states two six-rings per molecule.
    for tag, sizes in (("one", sizes1), ("two", sizes2)):
        version = "dseams1" if tag == "one" else "dseams2"
        for n in sizes:
            six = int(read_counts(DATA / f"{version}-ic-{n}.err")["six_rings"])
            if six != SIX_RINGS_PER_MOLECULE * n:
                raise SystemExit(f"{version} finds {six} six-rings on {n} molecules")
    for n in sizes2:
        b = [r[PATH] for r in samples[("two", n)]]
        if n in sizes1:
            a = [r[PATH] for r in samples[("one", n)]]
            rlo, rhi = bounds[(n, "path")]
            left = f"{len(a)} & {cells(a)}"
            ratio = (f"\\num{{{sig(median(a) / median(b))}}} & "
                     f"[\\num{{{sig(rlo)}}}, \\num{{{sig(rhi)}}}]")
        else:
            left, ratio = not_run(SAMPLE_COLUMNS), not_run(RATIO_COLUMNS)
        lines.append(f"\\num{{{n}}} & {left} & {len(b)} & "
                     f"{cells(b)} & {ratio}\\\\")
    mw = DATA / "dseams2-mw.txt"
    cols, rows = read_samples([mw])
    assert cols == ["ring_ms", "hc_ms", "ddc_ms", "path_ms"], cols
    for j, stage in enumerate(STAGES):
        define_median(f"v@mw@{stage}@two", [r[j] for r in rows])
    counts = read_counts(DATA / "dseams2-mw.err")
    for key, value in counts.items():
        define(f"v@mw@count@{key}", f"\\num{{{value}}}")
    b = [r[PATH] for r in rows]
    atoms = next(line.split()[1] for line in mw.read_text().splitlines()
                 if line.startswith("atoms "))
    if int(counts["six_rings"]) != SIX_RINGS_PER_MOLECULE * int(atoms):
        raise SystemExit(f"dseams2 finds {counts['six_rings']} six-rings on "
                         f"{atoms} mW molecules")
    define("v@mw@sites", f"\\num{{{atoms}}}")
    lines.append("\\addlinespace")
    lines.append(f"\\num{{{atoms}}}\\textsuperscript{{b}} & {not_run(SAMPLE_COLUMNS)} & "
                 f"{len(b)} & {cells(b)} & {not_run(RATIO_COLUMNS)}\\\\")
    write_table("tab-versions.tex", lines)
    top = common[-1]
    for stage in STAGES:
        for tag in ("one", "two"):
            xs = [r[STAGES.index(stage)] for r in samples[(tag, top)]]
            define_median(f"v@top@{stage}@{tag}", xs)
        lo, hi = bounds[(top, stage)]
        a = [r[STAGES.index(stage)] for r in samples[("one", top)]]
        b = [r[STAGES.index(stage)] for r in samples[("two", top)]]
        r = median(a) / median(b)
        define(f"v@top@{stage}@ratio", sig(r))
        define(f"v@top@{stage}@ratio@lo", sig(lo))
        define(f"v@top@{stage}@ratio@hi", sig(hi))
        define(f"v@top@{stage}@ratio@rel", relative(r, lo, hi))
        xs = [r[STAGES.index(stage)] for r in samples[("one", top)]]
        define(f"v@top@{stage}@one@s", sig(median(xs) / MS_PER_S))
    for tag, version in (("one", "dseams1"), ("two", "dseams2")):
        for key, value in read_counts(DATA / f"{version}-ic-{top}.err").items():
            define(f"v@top@count@{tag}@{key}", f"\\num{{{value}}}")
    rung = HEADLINE_RUNG
    one = median([r[PATH] for r in samples[("one", rung)]])
    two = median([r[PATH] for r in samples[("two", rung)]])
    define("hl@n", f"\\num{{{rung}}}")
    define("hl@one@s", sig(one / MS_PER_S, HEADLINE_FIGURES))
    define("hl@two@ms", sig(two, HEADLINE_FIGURES))
    define("hl@ratio", sig(one / two, HEADLINE_FIGURES))
    for j, stage in enumerate(STAGES):
        define(f"v@{rung}@{stage}@one@s",
               sig(median([r[j] for r in samples[("one", rung)]]) / MS_PER_S))
    define("hl@top@two@ms",
           sig(median([r[PATH] for r in samples[("two", top)]]), HEADLINE_FIGURES))
    names = {"ring": "Rings", "hc": "HC", "ddc": "DDC", "path": "Whole path"}
    n = common[-1]
    stage_lines = []
    for j, stage in enumerate(STAGES):
        a = [r[j] for r in samples[("one", n)]]
        b = [r[j] for r in samples[("two", n)]]
        rlo, rhi = bounds[(n, stage)]
        stage_lines.append(
            f"{names[stage]} & {cells(a)} & {cells(b)} & "
            f"\\num{{{sig(median(a) / median(b))}}} & [\\num{{{sig(rlo)}}}, "
            f"\\num{{{sig(rhi)}}}]\\\\")
    write_table("tab-stages.tex", stage_lines)
    for tag, sizes in (("one", sizes1), ("two", sizes2)):
        for j, stage in enumerate(STAGES):
            pts = [(n, median([r[j] for r in samples[(tag, n)]]))
                   for n in sizes]
            define(f"v@slope@{tag}@{stage}", sig(slope(pts), 2))
            out = ["sites median lo hi q1 q3"]
            for n in sizes:
                xs = [r[j] for r in samples[(tag, n)]]
                out.append(plot_row(n, xs))
            (DATA / f"plot-versions-{tag}-{stage}.dat").write_text(
                "\n".join(out) + "\n")
            plot_calls[f"plot-versions-{tag}-{stage}.dat"] = [
                len(samples[(tag, n)]) for n in sizes]
    # The main text gives the stages of the headline rung, one column each.
    head = []
    for tag, label in (("one", "d-SEAMS~1"), ("two", "d-SEAMS")):
        stages = [[r[j] for r in samples[(tag, rung)]] for j in range(len(STAGES))]
        found = [(median(xs), *median_interval(xs), *quartiles(xs)) for xs in stages]
        head.append(f"{label} & median & " + " & ".join(
            f"\\num{{{sig(m)}}}" for m, *_ in found) + "\\\\")
        head.append(f" & {100 * COVERAGE:g}\\% interval & " + " & ".join(
            f"[\\num{{{sig(lo)}}}, \\num{{{sig(hi)}}}]" for _, lo, hi, *_ in found)
            + "\\\\")
        head.append(" & relative error & " + " & ".join(
            f"\\qty{{{relative(m, lo, hi)}}}{{\\percent}}"
            for m, lo, hi, *_ in found) + "\\\\")
        head.append(" & quartiles & " + " & ".join(
            span(q1, q3) for *_, q1, q3 in found) + "\\\\")
        head.append("\\addlinespace")
    ratios = [median([r[j] for r in samples[("one", rung)]]) /
              median([r[j] for r in samples[("two", rung)]])
              for j in range(len(STAGES))]
    head.append("Ratio & median & " + " & ".join(
        f"\\num{{{sig(r)}}}" for r in ratios) + "\\\\")
    head.append(f" & {100 * COVERAGE:g}\\% interval & " + " & ".join(
        f"[\\num{{{sig(bounds[(rung, s)][0])}}}, \\num{{{sig(bounds[(rung, s)][1])}}}]"
        for s in STAGES) + "\\\\")
    head.append("\\addlinespace")
    for tag, label in (("one", "d-SEAMS~1"), ("two", "d-SEAMS")):
        sizes = sizes1 if tag == "one" else sizes2
        fits = [slope([(n, median([r[j] for r in samples[(tag, n)]]))
                       for n in sizes]) for j in range(len(STAGES))]
        first = "Fitted exponent" if tag == "one" else ""
        head.append(f"{first} & {label} & " + " & ".join(
            f"\\num{{{sig(f, 2)}}}" for f in fits) + "\\\\")
    write_table("tab-cost.tex", head)


CHURN_TIMES = ("ringInc_ms", "ringBatch_ms", "affilInc_ms", "affilBatch_ms")


def churn():
    """The update sweep; the rows print medians, the calls lines every call."""
    lines = (DATA / "churn-sweep.txt").read_text().splitlines()
    header = next(l for l in lines if l.startswith("# radius"))
    names = header[2:].split()
    rows = [dict(zip(names, l.split())) for l in lines
            if l and not l.startswith(("#", "calls "))]
    calls = {}
    for line in lines:
        if line.startswith("calls "):
            _, radius, key, *values = line.split()
            calls[(radius, key)] = [float(v) for v in values]
    meta = lines[0].split()
    reps = int(meta[meta.index("reps") + 1])
    nop = int(meta[meta.index("nop") + 1])
    define("c@reps", str(reps))
    define("c@nop", str(nop))
    for row in rows:
        for key in CHURN_TIMES:
            xs = calls[(row["radius"], key)]
            if len(xs) != reps or abs(median(xs) - float(row[key])) > ROW_TOLERANCE_MS:
                raise SystemExit(f"churn calls disagree with the row: {row['radius']} {key}")
            row[key] = median(xs)
            row[key + "@q"] = quartiles(xs)
    plot = ["churn ring_batch ring_inc affil_batch affil_inc ring_batch_q1 "
            "ring_batch_q3 ring_inc_q1 ring_inc_q3 affil_batch_q1 affil_batch_q3 "
            "affil_inc_q1 affil_inc_q3"]
    table = []
    dag = "\\textsuperscript{\\dag}"
    for row in rows:
        if int(row["bondChurn"]) == 0:
            continue
        order = ("ringBatch_ms", "ringInc_ms", "affilBatch_ms", "affilInc_ms")
        plot.append(" ".join([row["bondChurn"]] + [str(row[k]) for k in order] +
                             [str(q) for k in order for q in row[k + "@q"]]))
        marks = {"ringInc_ms": dag if int(row["ringRecomp"]) == nop else "",
                 "affilInc_ms": dag if row["affilReclass"] == row["nSix"] else ""}
        timing = " & ".join(f"\\num{{{sig(row[k])}}}{marks.get(k, '')} & "
                            f"{span(*row[k + '@q'])}" for k in CHURN_TIMES)
        label = f"\\num{{{row['radius']}}}" if row["radius"] != "-1" else "all"
        table.append(
            f"{label} & "
            f"\\num{{{row['changedAtoms']}}} & \\num{{{row['bondChurn']}}} & "
            f"{timing} & {row['ringsEqual']}/{row['affilEqual']}\\\\")
    (DATA / "plot-churn.dat").write_text("\n".join(plot) + "\n")
    write_table("tab-churn.tex", table)
    allequal = all(r["ringsEqual"] == "yes" and r["affilEqual"] == "yes"
                   for r in rows)
    define("c@allequal", "yes" if allequal else "no")
    moving = [r for r in rows if int(r["bondChurn"]) > 0]
    small = moving[0]
    define("c@small@churn", small["bondChurn"])
    define("c@small@atoms", small["changedAtoms"])
    for key in CHURN_TIMES:
        define(f"c@small@{key}", sig(float(small[key])))
        define(f"c@small@{key}@q1", sig(small[key + "@q"][0]))
        define(f"c@small@{key}@q3", sig(small[key + "@q"][1]))
    define("c@small@ringgain", sig(float(small["ringBatch_ms"]) /
                                   float(small["ringInc_ms"])))
    define("c@small@affgain", sig(float(small["affilBatch_ms"]) /
                                  float(small["affilInc_ms"])))
    define("hl@ringgain", sig(float(small["ringBatch_ms"]) /
                              float(small["ringInc_ms"]), 2))
    nop = int(meta[meta.index("nop") + 1])
    spliced = [r for r in moving if int(r["ringRecomp"]) < nop]
    reclassified = [r for r in moving
                    if int(r["affilReclass"]) < int(r["nSix"])]
    last = [r for r in spliced
            if float(r["ringInc_ms"]) < float(r["ringBatch_ms"])][-1]
    define("c@ringlast@churn", last["bondChurn"])
    define("c@ringlast@atoms", last["changedAtoms"])
    define("c@ringlast@inc", sig(float(last["ringInc_ms"])))
    define("c@ringlast@batch", sig(float(last["ringBatch_ms"])))
    alast = [r for r in reclassified
             if float(r["affilInc_ms"]) < float(r["affilBatch_ms"])][-1]
    define("c@afflast@churn", alast["bondChurn"])
    define("c@afflast@atoms", alast["changedAtoms"])
    define("c@afflast@inc", sig(float(alast["affilInc_ms"])))
    define("c@afflast@batch", sig(float(alast["affilBatch_ms"])))
    define("c@spliced", str(len(spliced)))
    define("c@reclassified", str(len(reclassified)))
    rfall = next(r for r in moving if int(r["ringRecomp"]) == nop)
    define("c@ringfall@churn", rfall["bondChurn"])
    define("c@ringfall@atoms", rfall["changedAtoms"])
    define("c@ringfall@perbond", sig(int(rfall["bondChurn"]) / nop, 2))
    afall = next(r for r in moving
                 if int(r["affilReclass"]) == int(r["nSix"]))
    define("c@afffall@churn", afall["bondChurn"])
    define("c@afffall@atoms", afall["changedAtoms"])
    fallbacks = ([(r["ringInc_ms"], r["ringBatch_ms"]) for r in moving
                  if int(r["ringRecomp"]) == nop] +
                 [(r["affilInc_ms"], r["affilBatch_ms"]) for r in moving
                  if r["affilReclass"] == r["nSix"]])
    # The main text gives the smallest change and the limits of both updates.
    update = []
    for name, inc, full, gain, fall in (
            ("rings", "ringInc_ms", "ringBatch_ms", last, rfall),
            ("cage labels", "affilInc_ms", "affilBatch_ms", alast, afall)):
        update.append(
            f"{name} & \\num{{{sig(small[full])}}} & {span(*small[full + '@q'])} & "
            f"\\num{{{sig(small[inc])}}} & {span(*small[inc + '@q'])} & "
            f"\\num{{{sig(small[full] / small[inc])}}} & "
            f"\\num{{{gain['bondChurn']}}} & \\num{{{fall['bondChurn']}}}\\\\")
    write_table("tab-update.tex", update)
    define("c@fallbacks", str(len(fallbacks)))
    define("c@fallback@below",
           "yes" if all(inc < full for inc, full in fallbacks) else "no")
    define("c@fallback@worst", sig(max(inc / full for inc, full in fallbacks)))
    full = rows[-1]
    define("c@full@churn", full["bondChurn"])
    define("c@full@atoms", full["changedAtoms"])
    define("c@nsix", rows[0]["nSix"])


H2H_SYSTEMS = ("ic", "ih", "liquid", "confined")
SYSTEM_NAMES = {"ic": "Ice~Ic", "ih": "Ice~Ih", "liquid": "Liquid",
                "confined": "Confined"}


def h2h():
    timings = list(csv.DictReader((H2H / "timings.tsv").open(), delimiter="\t"))
    comps = list(csv.DictReader((H2H / "comparisons.tsv").open(), delimiter="\t"))
    agree = list(csv.DictReader((H2H / "agreements.tsv").open(), delimiter="\t"))
    census = list(csv.DictReader((H2H / "ring-census.tsv").open(), delimiter="\t"))
    errors = [json.loads(l) for l in (H2H / "errors.jsonl").read_text().splitlines()
              if l.strip()]
    timing = {}
    for row in timings:
        key = (row["system"], row["scale"], row["method"])
        timing[key] = row
        base = f"h@{row['system']}@{row['scale']}@{row['method']}"
        m = float(row["median_ms"])
        lo, hi = float(row["ci95_low_ms"]), float(row["ci95_high_ms"])
        define(base, sig(m))
        define(base + "@lo", sig(lo))
        define(base + "@hi", sig(hi))
        define(base + "@rel", relative(m, lo, hi))
        define(base + "@q1", sig(float(row["q1_ms"])))
        define(base + "@q3", sig(float(row["q3_ms"])))
        define(base + "@n", row["repeats"])
        define(f"h@{row['system']}@{row['scale']}@sites", f"\\num{{{row['molecules']}}}")
    ratio = {}
    for row in comps:
        key = (row["system"], row["scale"], row["other"])
        ratio[key] = row
        base = f"r@{row['system']}@{row['scale']}@{row['other']}"
        r = float(row["other_over_dseams"])
        lo, hi = float(row["ratio_ci95_low"]), float(row["ratio_ci95_high"])
        define(base, sig(r))
        define(base + "@lo", sig(lo))
        define(base + "@hi", sig(hi))
        define(base + "@rel", relative(r, lo, hi))
        define(base + "@p", sig(float(row["p_holm"]), 2))
    cycless = {sig(float(ratio[(s, "1e3", "genice-rings")]["other_over_dseams"]), 2)
               for s in ("ic", "ih")}
    if len(cycless) != 1:
        raise SystemExit(f"cycless ratios on the crystals round apart: {cycless}")
    define("hl@cycless", cycless.pop())
    supported = [r for r in comps if r["equivalent_workload"] == "True"]
    define("h@equivalent@maxratio",
           sig(max(float(r["other_over_dseams"]) for r in supported)))
    define("h@nonequivalent", str(len(comps) - len(supported)))
    define("h@comparisons", str(len(comps)))
    best = {}
    for r in supported:
        if r["other"] in ("ovito-chill", "palmer-chill", "freud-q6", "pyscal-q6",
                          "matscipy-rings-kernel", "matscipy-rings-public"):
            key = r["other"]
            value = float(r["other_over_dseams"])
            lo, hi = best.get(key, (value, value))
            best[key] = (min(lo, value), max(hi, value))
    for key, (lo, hi) in best.items():
        define(f"h@range@{key}@min", sig(lo))
        define(f"h@range@{key}@max", sig(hi))
    matscipy = [v for k, v in best.items() if k.startswith("matscipy-rings-")]
    if matscipy:
        define("h@range@matscipy@min", sig(min(lo for lo, _ in matscipy)))
        define("h@range@matscipy@max", sig(max(hi for _, hi in matscipy)))
    local = [r for r in comps if r["operation"] == "ice-polymorph"]
    for other in ("ovito-ptm", "ovito-diamond", "ovito-cna", "pyscal-diamond"):
        values = [float(r["other_over_dseams"]) for r in local
                  if r["other"] == other]
        define(f"h@local@{other}@min", sig(min(values)))
        define(f"h@local@{other}@max", sig(max(values)))
        define(f"h@local@{other}@invmin", sig(1.0 / max(values)))
        define(f"h@local@{other}@invmax", sig(1.0 / min(values)))
    define("hl@local@max",
           sig(max(1.0 / float(r["other_over_dseams"]) for r in local), 2))
    define("h@equivalent", str(len(supported)))
    define("h@equivalent@ordered",
           str(sum(1 for r in supported if r["ordering_supported"] == "True")))
    define("h@equivalent@minratio",
           sig(min(float(r["other_over_dseams"]) for r in supported)))
    define("h@equivalent@maxp",
           sig(max(float(r["p_holm"]) for r in supported), 2))
    define("h@records", str(len(timings)))
    define("h@errors", str(len(errors)))
    for row in agree:
        base = f"a@{row['system']}@{row['scale']}@{row['left']}@{row['right']}"
        define(base, f"\\num{{{row['agreed']}}}")
        total = int(row["agreed"]) + int(row["disagreed"])
        define(base + "@total", f"\\num{{{total}}}")
        define(base + "@pct", sig(100.0 * int(row["agreed"]) / total, 3))
        define(base + "@off", sig(100.0 * int(row["disagreed"]) / total, 2))
    for row in census:
        define(f"k@{row['system']}@{row['scale']}@{row['method']}",
               ", ".join(f"\\num{{{row[f'rings_{s}']}}}" for s in RING_SIZES))
        define(f"k@{row['system']}@{row['scale']}@{row['method']}@six",
               f"\\num{{{row['rings_6']}}}")
    by_method = {}
    for err in errors:
        by_method.setdefault(err["method"], []).append(err)
    for method, errs in by_method.items():
        define(f"e@{method}", str(len(errs)))
    plots = {
        "chill": ("ic", ["dseams-chill", "ovito-chill", "palmer-chill"]),
        "q6": ("ic", ["dseams-q6", "freud-q6", "pyscal-q6"]),
        "rings": ("liquid", ["dseams-rings", "matscipy-rings-kernel"]),
        "seeded": ("ic", ["dseams-seeded", "ovito-ptm", "ovito-diamond",
                          "ovito-cna", "pyscal-diamond"]),
    }
    for name, (system, methods) in plots.items():
        for method in methods:
            rows = [timing[k] for k in sorted(timing, key=lambda k: float(k[1]))
                    if k[0] == system and k[2] == method]
            out = ["sites median lo hi q1 q3"] + [
                f"{r['molecules']} {r['median_ms']} {r['ci95_low_ms']} "
                f"{r['ci95_high_ms']} {r['q1_ms']} {r['q3_ms']}" for r in rows]
            (DATA / f"plot-{name}-{method}.dat").write_text("\n".join(out) + "\n")
            plot_calls[f"plot-{name}-{method}.dat"] = [int(r["repeats"]) for r in rows]
    # The method identifier names the operation, so a frame is one group of
    # rows in the order of the operations below.
    ladder = []
    operations = ("chill", "ice-liquid", "rings", "ice-polymorph", "rings-process")
    for system in H2H_SYSTEMS:
        for scale in ("1e3", "1e4", "1e5", "1e6"):
            frame = [r for op in operations for r in timings
                     if r["system"] == system and r["scale"] == scale and
                     r["operation"] == op]
            if not frame:
                continue
            if ladder:
                ladder.append("\\addlinespace")
            for i, r in enumerate(frame):
                comp = ratio.get((system, scale, r["method"]))
                rtext = (f"\\num{{{sig(float(comp['other_over_dseams']))}}}"
                         if comp else "")
                head = (f"{SYSTEM_NAMES[system]} & \\num{{{r['molecules']}}}"
                        if i == 0 else " & ")
                ladder.append(
                    f"{head} & \\texttt{{{r['method']}}} & {r['repeats']} & "
                    f"\\num{{{sig(float(r['median_ms']))}}} & "
                    f"[\\num{{{sig(float(r['ci95_low_ms']))}}}, "
                    f"\\num{{{sig(float(r['ci95_high_ms']))}}}] & "
                    f"{span(float(r['q1_ms']), float(r['q3_ms']))} & {rtext}\\\\")
    write_table("tab-ladder.tex", ladder)
    cens = []
    for err in errors:
        cens.append((err["system"], err["scale"], err["method"],
                     err.get("status", ""), err.get("error", "")))
    summary = {}
    for system, scale, method, status, text in cens:
        kind = "timeout" if "timeout" in (status + text).lower() or \
            "timed out" in text.lower() else (
                "signal 11" if "-11" in text else "error")
        name = SYSTEM_NAMES[system]
        summary.setdefault((method, kind), []).append(
            f"{name[0].lower()}{name[1:]} \\(10^{{{scale[2:]}}}\\)")
    lines = []
    for (method, kind), where in sorted(summary.items()):
        lines.append(f"\\texttt{{{method}}} & {kind} & {len(where)} & "
                     f"{', '.join(where)}\\\\")
    write_table("tab-censored.tex", lines)


def worked():
    rows = list(csv.DictReader((DATA / "worked_example.tsv").open(), delimiter="\t"))
    for row in rows:
        name = row["frame"]
        define(f"w@{name}@sites", f"\\num{{{row['sites']}}}")
        rings = row["rings_3_to_6"].split(",")
        for size, count in zip(RING_SIZES, rings):
            define(f"w@{name}@ring{size}", f"\\num{{{count}}}")
        hexa, cubic, water, both = row["cage_hex_cubic_water_both"].split(",")
        define(f"w@{name}@ih", f"\\num{{{hexa}}}")
        define(f"w@{name}@ic", f"\\num{{{cubic}}}")
        define(f"w@{name}@water", f"\\num{{{water}}}")
        prisms = row["prisms_4_5_6"].split(",")
        heights = row["height_percent_4_5_6"].split(",")
        for size, count, height in zip(PRISM_SIZES, prisms, heights):
            define(f"w@{name}@prism{size}", f"\\num{{{count}}}")
            define(f"w@{name}@height{size}", sig(float(height), 4))
    sized = {row["frame"] for row in rows}
    ideal = DATA / "ideal-cages.txt"
    names = {"4:6,6:2": "hc", "6:7": "ddc", "5:12": "512",
             "5:12,6:2": "51262", "5:12,6:4": "51264"}
    frame = None
    for line in ideal.read_text().splitlines():
        parts = line.split()
        if parts and parts[0] == "#":
            stem = next(p for p in parts if p.endswith(".lammpstrj"))
            frame = pathlib.Path(stem).stem.removeprefix("genice_")
        elif parts and parts[0] == "nop":
            fields = dict(zip(parts[0::2], parts[1::2]))
            if frame not in sized:
                sized.add(frame)
                define(f"w@{frame}@sites", f"\\num{{{fields['nop']}}}")
            label = names[fields["signature"]]
            define(f"w@{frame}@cages@{label}", f"\\num{{{fields['cages']}}}")
            define(f"w@{frame}@cageatoms@{label}", f"\\num{{{fields['atoms']}}}")
    update = DATA / "worked_update.tsv"
    for row in csv.DictReader(update.open(), delimiter="\t"):
        value = row["value"]
        if "," in value:
            for size, count in zip(RING_SIZES, value.split(",")):
                define(f"u@{row['quantity']}@{size}", f"\\num{{{count}}}")
        elif value.isdigit():
            define(f"u@{row['quantity']}", f"\\num{{{value}}}")


def nucleation():
    """Cage labels along the nucleation trajectory, carried and rebuilt."""
    labels = {}
    for graph in ("cutoff", "seeded"):
        rows = []
        for line in (DATA / f"nucleation-{graph}-carry.txt").read_text().splitlines():
            parts = line.split()
            if parts and parts[0] == "nop":
                fields = dict(zip(parts[0::2], parts[1::2]))
                rows.append(fields)
        labels[graph] = rows
    frames = len(labels["cutoff"])
    define("n@frames", str(frames))
    define("n@sites", f"\\num{{{labels['cutoff'][0]['nop']}}}")
    out = ["frame cut_hex cut_cubic cut_both seed_hex seed_cubic seed_both"]
    for i, (a, b) in enumerate(zip(labels["cutoff"], labels["seeded"]), 1):
        out.append(" ".join([str(i), a["hexagonal"], a["cubic"], a["both"],
                             b["hexagonal"], b["cubic"], b["both"]]))
    (DATA / "plot-nucleation.dat").write_text("\n".join(out) + "\n")
    for graph, rows in labels.items():
        for tag, row in (("first", rows[0]), ("last", rows[-1])):
            for key in ("hexagonal", "cubic", "both", "water"):
                define(f"n@{graph}@{tag}@{key}", f"\\num{{{row[key]}}}")
            ice = sum(int(row[k]) for k in ("hexagonal", "cubic", "both"))
            define(f"n@{graph}@{tag}@ice", f"\\num{{{ice}}}")
            cubhex = int(row["cubic"]) / int(row["hexagonal"])
            define(f"n@{graph}@{tag}@cubhex",
                   sig(cubhex, 3 if f"{cubhex:e}".startswith("1") else 2))
            define(f"n@{graph}@{tag}@both@pct", str(round(100.0 * int(row["both"]) / ice)))
        totals = [sum(int(r[k]) for k in ("hexagonal", "cubic", "both")) for r in rows]
        top = max(range(len(totals)), key=totals.__getitem__)
        define(f"n@{graph}@max@ice", f"\\num{{{totals[top]}}}")
        define(f"n@{graph}@max@frame", str(top + 1))
        gains = {k: int(rows[-1][k]) - int(rows[0][k])
                 for k in ("hexagonal", "cubic", "both")}
        gained = sum(gains.values())
        define(f"n@{graph}@gain@ice", f"\\num{{{gained}}}")
        for key, value in gains.items():
            define(f"n@{graph}@gain@{key}", f"\\num{{{value}}}")
            define(f"n@{graph}@gain@{key}@pct", str(round(100.0 * value / gained)))
        define(f"n@{graph}@gain@cubicboth@pct",
               str(round(100.0 * (gains["cubic"] + gains["both"]) / gained)))
    table = []
    classes = ("hexagonal", "cubic", "both")
    for graph, name in (("cutoff", "Cutoff"), ("seeded", "Seeded")):
        rows = labels[graph]
        totals = [sum(int(r[k]) for k in classes) for r in rows]
        top = max(range(len(totals)), key=totals.__getitem__)
        sets = [("initial", {k: int(rows[0][k]) for k in classes}),
                ("last", {k: int(rows[-1][k]) for k in classes}),
                (f"largest, frame {top + 1}", {k: int(rows[top][k]) for k in classes}),
                ("gain", {k: int(rows[-1][k]) - int(rows[0][k]) for k in classes})]
        for i, (label, c) in enumerate(sets):
            ice = sum(c.values())
            table.append(
                f"{name if i == 0 else ''} & {label} & "
                + " & ".join(f"\\num{{{c[k]}}}" for k in classes)
                + f" & \\num{{{ice}}} & \\num{{{sig(c['cubic'] / c['hexagonal'], 3)}}}"
                f" & \\qty{{{round(100.0 * c['both'] / ice)}}}{{\\percent}}\\\\")
        table.append("\\addlinespace")
    write_table("tab-nucleation-cages.tex", table[:-1])
    shares = {g: {k: 100.0 * (int(labels[g][-1][k]) - int(labels[g][0][k])) /
                  sum(int(labels[g][-1][j]) - int(labels[g][0][j]) for j in classes)
                  for k in classes} for g in labels}
    define("n@gain@maxdiff@pct",
           sig(max(abs(shares["cutoff"][k] - shares["seeded"][k]) for k in classes), 2))
    host = (DATA / "nucleation-host.txt").read_text().split("\n")
    for line in host:
        parts = line.split()
        if len(parts) == 3 and parts[1] == "identical":
            define(f"n@{parts[0]}@identical", parts[2])
    times = {}
    for line in (DATA / "nucleation-timings.txt").read_text().splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0] in ("cutoff", "seeded"):
            times.setdefault((parts[0], parts[1]), []).append(float(parts[3]))
    rng = random.Random(SEED)
    for (graph, mode), xs in times.items():
        define_median(f"n@{graph}@{mode}", xs)
        define(f"n@{graph}@{mode}@frame", sig(median(xs) / frames))
        q1, q3 = quartiles(xs)
        define(f"n@{graph}@{mode}@frame@q1", sig(q1 / frames))
        define(f"n@{graph}@{mode}@frame@q3", sig(q3 / frames))
    for graph in ("cutoff", "seeded"):
        a = times[(graph, "rebuild")]
        b = times[(graph, "carry")]
        lo, hi = ratio_interval(a, b, rng)
        r = median(a) / median(b)
        define(f"n@{graph}@ratio", sig(r))
        define(f"n@{graph}@ratio@lo", sig(lo))
        define(f"n@{graph}@ratio@hi", sig(hi))


def neighbours():
    path = DATA / "h2h" / "chill-neighbours.txt"
    for line in path.read_text().splitlines():
        key, value = line.split()
        define(f"x@{key}", f"\\num{{{value}}}")


def nucleation_churn():
    """Changed bonds and atoms between saved nucleation frames."""
    fields = dict(line.split() for line in
                  (DATA / "nucleation-churn.txt").read_text().splitlines())
    for key in ("transitions", "molecules", "over_atom_rule",
                "over_bond_rule", "over_either_rule", "median_changed_atoms"):
        define(f"m@{key}", f"\\num{{{fields[key]}}}")
    for key in ("median", "min", "max", "q1", "q3"):
        define(f"m@{key}@bonds", sig(float(fields[f"{key}_bonds_per_molecule"])))
    for key in ("q1", "q3"):
        define(f"m@{key}@atoms", sig(float(fields[f"{key}_changed_atoms"])))
    define("m@median@atoms@pct", sig(100.0 * int(fields["median_changed_atoms"]) /
                                     int(fields["molecules"]), 2))


def ovito_rings():
    """d-SEAMS rings against OVITO FindRings, paired in one session."""
    rows = list(csv.DictReader((H2H / "ovito-rings.tsv").open(), delimiter="\t"))
    rng = random.Random(SEED)
    samples = {}
    census = {}
    for row in rows:
        system, scale, method = row["system"], row["scale"], row["method"]
        xs = [float(x) for x in row["samples_ms"].split(",")]
        samples[(system, scale, method)] = xs
        census.setdefault((system, scale), set()).add(row["rings_3_to_6"])
        define(f"o@{system}@{scale}@sites", f"\\num{{{row['molecules']}}}")
        define_median(f"o@{system}@{scale}@{method}", xs)
    ratios = []
    separated = True
    for system in H2H_SYSTEMS:
        out = {m: ["sites median lo hi q1 q3"] for m in ("dseams-rings", "ovito-rings")}
        for scale in ("1e3", "1e4", "1e5", "1e6"):
            a = samples[(system, scale, "ovito-rings")]
            b = samples[(system, scale, "dseams-rings")]
            lo, hi = ratio_interval(a, b, rng)
            r = median(a) / median(b)
            base = f"o@{system}@{scale}@ratio"
            define(base, sig(r))
            define(base + "@lo", sig(lo))
            define(base + "@hi", sig(hi))
            define(base + "@rel", relative(r, lo, hi))
            ratios.append(r)
            separated = separated and min(a) > max(b)
            sites = next(x["molecules"] for x in rows
                         if x["system"] == system and x["scale"] == scale)
            for method, xs in (("dseams-rings", b), ("ovito-rings", a)):
                out[method].append(plot_row(sites, xs))
                plot_calls.setdefault(f"plot-ovrings-{system}-{method}.dat",
                                      []).append(len(xs))
        for method, lines in out.items():
            (DATA / f"plot-ovrings-{system}-{method}.dat").write_text(
                "\n".join(lines) + "\n")
    define("o@ratio@min", sig(min(ratios)))
    define("o@ratio@max", sig(max(ratios)))
    define("o@separated", "yes" if separated else "no")
    define("o@frames", str(len(census)))
    define("o@comparisons", str(len(ratios)))
    define("o@census@equal",
           "yes" if all(len(v) == 1 for v in census.values()) else "no")
    return samples


H2H_ROWS = (
    ("CHILL+", "dseams-chill", (("OVITO CHILL+", "ovito-chill"),
                                ("CHILL+ library", "palmer-chill"))),
    ("\\(q_6\\)", "dseams-q6", (("freud", "freud-q6"), ("pyscal", "pyscal-q6"))),
    ("cage labels", "dseams-seeded", (("OVITO PTM", "ovito-ptm"),
                                    ("OVITO diamond", "ovito-diamond"),
                                    ("OVITO CNA", "ovito-cna"),
                                    ("pyscal diamond", "pyscal-diamond"))),
)


def h2h_table(rings, scale="1e5"):
    """d-SEAMS medians and comparator ratios at one scale."""
    timings = list(csv.DictReader((H2H / "timings.tsv").open(), delimiter="\t"))
    comps = list(csv.DictReader((H2H / "comparisons.tsv").open(), delimiter="\t"))
    errors = [json.loads(l) for l in (H2H / "errors.jsonl").read_text().splitlines()
              if l.strip()]
    time = {(r["system"], r["method"]): r for r in timings if r["scale"] == scale}
    ratio = {(r["system"], r["other"]): r for r in comps if r["scale"] == scale}
    failed = {(e["system"], e["method"]): e for e in errors if e["scale"] == scale}
    rels = []

    def failure(system, method):
        err = failed.get((system, method))
        if err is None:
            raise SystemExit(f"no timing and no failure: {system} {scale} {method}")
        text = (err.get("status", "") + " " + err.get("error", "")).lower()
        if "timeout" in text or "timed out" in text:
            return "{timeout\\textsuperscript{c}}"
        if "-11" in text:
            return "{signal 11\\textsuperscript{a}}"
        return "{error\\textsuperscript{b}}"

    def quartile_row(spans):
        return " & quartiles & " + " & ".join(spans) + "\\\\"

    lines = []
    for label, ours, others in H2H_ROWS:
        cells = [f"\\num{{{sig(float(time[(s, ours)]['median_ms']))}}}"
                 for s in H2H_SYSTEMS]
        lines.append(f"{label} & d-SEAMS (\\unit{{\\ms}}) & " + " & ".join(cells) + "\\\\")
        lines.append(quartile_row(
            span(float(time[(s, ours)]["q1_ms"]), float(time[(s, ours)]["q3_ms"]))
            for s in H2H_SYSTEMS))
        for name, other in others:
            cells = []
            for s in H2H_SYSTEMS:
                row = ratio.get((s, other))
                if row is None:
                    cells.append(failure(s, other))
                    continue
                r = float(row["other_over_dseams"])
                lo, hi = float(row["ratio_ci95_low"]), float(row["ratio_ci95_high"])
                rels.append(float(relative(r, lo, hi)))
                cells.append(f"\\num{{{sig(r)}}}")
            lines.append(f" & {name} & " + " & ".join(cells) + "\\\\")
        lines.append("\\addlinespace")
    rng = random.Random(SEED)
    cells = [f"\\num{{{sig(median(rings[(s, scale, 'dseams-rings')]))}}}"
             for s in H2H_SYSTEMS]
    lines.append("rings & d-SEAMS (\\unit{\\ms}) & " + " & ".join(cells) + "\\\\")
    lines.append(quartile_row(span(*quartiles(rings[(s, scale, "dseams-rings")]))
                              for s in H2H_SYSTEMS))
    cells = []
    for s in H2H_SYSTEMS:
        a = rings[(s, scale, "ovito-rings")]
        b = rings[(s, scale, "dseams-rings")]
        lo, hi = ratio_interval(a, b, rng)
        r = median(a) / median(b)
        rels.append(float(relative(r, lo, hi)))
        cells.append(f"\\num{{{sig(r)}}}")
    lines.append(" & OVITO rings & " + " & ".join(cells) + "\\\\")
    for name, method in (("matscipy, public", "matscipy-rings-public"),
                         ("matscipy, compiled", "matscipy-rings-kernel"),
                         ("cycless", "genice-rings"), ("RINGS", "rings-1.3.5")):
        if any((s, method) in ratio for s in H2H_SYSTEMS):
            raise SystemExit(f"{method} completed at {scale}; the table lists only failures")
        cells = [failure(s, method) for s in H2H_SYSTEMS]
        lines.append(f" & {name} & " + " & ".join(cells) + "\\\\")
    write_table("tab-h2h.tex", lines)
    define("t@maxrel", sig(max(rels), 2))


MATCHED = (("OVITO", "CHILL+", "ovito-chill"),
           ("CHILL+ library", "CHILL+", "palmer-chill"),
           ("freud", "\\(q_6\\)", "freud-q6"),
           ("pyscal", "\\(q_6\\)", "pyscal-q6"),
           ("matscipy, public", "primitive rings", "matscipy-rings-public"),
           ("matscipy, compiled", "primitive rings", "matscipy-rings-kernel"))
APART = (("OVITO PTM", "cage labels", "ovito-ptm"),
         ("OVITO diamond", "cage labels", "ovito-diamond"),
         ("OVITO CNA", "cage labels", "ovito-cna"),
         ("pyscal diamond", "cage labels", "pyscal-diamond"),
         ("cycless", "rings", "genice-rings"))


def pvalue(p):
    return sig(p, 2) if p >= PVALUE_FIXED else f"{p:.1e}"


def ranges(rings):
    """Smallest and largest ratio of each code to d-SEAMS, for the main text."""
    comps = list(csv.DictReader((H2H / "comparisons.tsv").open(), delimiter="\t"))
    supported = {r["other"] for r in comps if r["equivalent_workload"] == "True"}
    if supported != {m for *_, m in MATCHED}:
        raise SystemExit(f"matched codes differ from the records: {sorted(supported)}")

    def line(code, operation, values, p):
        return (f"{code} & {operation} & {len(values)} & \\num{{{sig(min(values))}}} & "
                f"\\num{{{sig(max(values))}}} & {p}\\\\")

    def campaign(code, operation, method):
        chosen = [r for r in comps if r["other"] == method]
        values = [float(r["other_over_dseams"]) for r in chosen]
        p = max(float(r["p_holm"]) for r in chosen)
        return values, line(code, operation, values, f"\\num{{{pvalue(p)}}}")

    lines = ["\\multicolumn{6}{@{}l}{Matched pairs}\\\\"]
    matched = []
    for code, operation, method in MATCHED:
        values, text = campaign(code, operation, method)
        matched.extend(values)
        lines.append(text)
    scales = sorted({k[1] for k in rings}, key=float)
    ovito = [median(rings[(s, c, "ovito-rings")]) / median(rings[(s, c, "dseams-rings")])
             for s in H2H_SYSTEMS for c in scales]
    matched.extend(ovito)
    lines.append(line("OVITO", "primitive rings", ovito, "\\textsuperscript{a}"))
    lines.append("\\addlinespace")
    lines.append("\\multicolumn{6}{@{}l}{Reported apart}\\\\")
    faster = 0
    for code, operation, method in APART:
        values, text = campaign(code, operation, method)
        lines.append(text)
        faster += method != "genice-rings" and max(values) < 1.0
    write_table("tab-ranges.tex", lines)
    if min(matched) <= 1.0:
        raise SystemExit("a matched code is not slower than d-SEAMS")
    define("r@matched", str(len(matched)))
    define("r@local@faster", str(faster))


THRESHOLDS = (
    ("ringatoms", "franzblau.cpp", r"double fullPassFraction = ([0-9.]+);"),
    ("ringbonds", "franzblau.cpp", r"double bondCrossover = ([0-9.]+);"),
    ("cagerows", "cage_affiliation.cpp", r"double batchFraction = ([0-9.]+);"),
    ("packcap", "cage_pack.cpp", r"static constexpr int kCap = ([0-9]+);"),
    ("edgelo", "include/internal/cage_enum.hpp",
     r"struct EdgeWindow \{[^}]*double lo = ([0-9.]+);"),
    ("edgehi", "include/internal/cage_enum.hpp",
     r"struct EdgeWindow \{[^}]*double hi = ([0-9.]+);"),
)


# Copies of the thresholds in the scripts, which must equal the defaults in src/.
THRESHOLD_COPIES = (
    ("ringatoms", "paper/tex/figures/nucleation_churn.py", r"^ATOM_FRACTION = ([0-9.]+)$"),
    ("ringbonds", "paper/tex/figures/nucleation_churn.py", r"^BOND_CROSSOVER = ([0-9.]+)$"),
)


def thresholds():
    """The fallback thresholds and fixed limits, read from their defaults in src/."""
    found_values = {}
    for name, path, pattern in THRESHOLDS:
        found = re.findall(pattern, (SOURCE / path).read_text())
        if len(found) != 1:
            raise SystemExit(f"{name}: {len(found)} matches in src/{path}")
        value = decimal.Decimal(found[0])
        found_values[name] = value
        define(f"thr@{name}", str(value))
        if value < 1:
            define(f"thr@{name}@pct", f"{(100 * value).normalize():f}")
    for name, path, pattern in THRESHOLD_COPIES:
        if read_setting(name, ((path, pattern),)) != found_values[name]:
            raise SystemExit(f"{path} copies {name} with another value than src/")


def read_setting(name, sources):
    """The one value that every source gives for a setting."""
    found = set()
    for path, pattern in sources:
        values = re.findall(pattern, (ROOT / path).read_text(), re.M)
        if not values:
            raise SystemExit(f"{name}: no match in {path}")
        found.update(decimal.Decimal(v) for v in values)
    if len(found) != 1:
        raise SystemExit(f"{name}: the files disagree: {sorted(found)}")
    return found.pop()


def hardware():
    """The machine records of the timings must name the same processor."""
    fields = []
    for path in (DATA / "host.txt", DATA / "churn-host.txt"):
        text = path.read_text()
        fields.append((re.search(r"^uname -m: (\S+)", text, re.M).group(1),
                       re.search(r"^nproc: (\d+)", text, re.M).group(1),
                       re.search(r"^model name\s*: (.+)$", text, re.M).group(1).strip(),
                       re.search(r"^cpu family\s*: (\d+)", text, re.M).group(1),
                       re.search(r"^model\s*: (\d+)", text, re.M).group(1)))
    lscpu = dict(line.split(":", 1) for line in
                 json.loads((H2H / "environment.json").read_text())["lscpu"].splitlines()
                 if ":" in line)
    lscpu = {k.strip(): v.strip() for k, v in lscpu.items()}
    fields.append((lscpu["Architecture"], lscpu["CPU(s)"], lscpu["Model name"],
                   lscpu["CPU family"], lscpu["Model"]))
    if len(set(fields)) != 1:
        raise SystemExit(f"the machine records differ: {fields}")
    # The text gives the architecture, the processor count and the model name;
    # the family and the model number only have to agree.
    arch, nproc, model = fields[0][:3]
    host = (DATA / "host.txt").read_text()
    define("hw@arch", arch.replace("_", "\\_"))
    define("hw@nproc", nproc)
    define("hw@model", model)
    define("hw@gcc", re.search(r"^c\+\+ \(GCC\) (\S+)", host, re.M).group(1))
    return fields[0]


def marks():
    """Largest spread of a plotted interval or quartile, as a factor of the median.

    The captions state that these spreads lie inside the marks of MARK_SPAN.
    """
    worst = {}
    for path in DATA.glob("plot-*.dat"):
        lines = path.read_text().split("\n")
        header = lines[0].split()
        family = path.stem.split("-")[1]
        for line in lines[1:]:
            if not line or header[0] not in ("sites", "churn"):
                continue
            row = dict(zip(header, map(float, line.split())))
            if header[0] == "sites":
                pairs = [("median", "lo", "hi"), ("median", "q1", "q3")]
            else:
                pairs = [(k, f"{k}_q1", f"{k}_q3") for k in
                         ("ring_batch", "ring_inc", "affil_batch", "affil_inc")]
            for m, lo, hi in pairs:
                spread = max(row[hi] / row[m], row[m] / row[lo])
                worst[family] = max(worst.get(family, 1.0), spread)
    for family, spread in worst.items():
        if spread > MARK_SPAN.get(family, MARK_SPAN["other"]):
            raise SystemExit(f"plot-{family}: a spread of {spread:.3f} leaves the marks")
    return worst


def figure_calls():
    """The fewest calls behind a mark of each figure of CALL_FIGURES."""
    for figure in CALL_FIGURES:
        tex = (pathlib.Path(__file__).parent / f"fig-{figure}.tex").read_text()
        plots = re.findall(r"\{figures/data/(plot-[^}]*\.dat)\}", tex)
        if not plots:
            raise SystemExit(f"fig-{figure}.tex plots no records")
        define(f"fig@{figure}@calls",
               word(min(calls for plot in plots for calls in plot_calls[plot])))


def main():
    versions()
    churn()
    h2h()
    worked()
    nucleation()
    neighbours()
    nucleation_churn()
    thresholds()
    rings = ovito_rings()
    h2h_table(rings)
    ranges(rings)
    machine = hardware()
    spread = marks()
    figure_calls()
    sys.stdout.write(f"machine {' | '.join(machine)}\n")
    sys.stdout.write("largest plotted spread " + ", ".join(
        f"{k} {v:.3f}" for k, v in sorted(spread.items())) + "\n")
    count = {}
    for line in macros:
        for key in ("h@equivalent", "o@comparisons"):
            if f"dsv@{key}\\endcsname" in line:
                count[key] = int(line.rsplit("{", 1)[1].rstrip("}"))
    define("hl@matched", str(count["h@equivalent"] + count["o@comparisons"]))
    body = {}
    for line in macros:
        name, _, rest = line.partition("dsv@")[2].partition("\\endcsname")
        if body.setdefault(name, rest[1:-1]) != rest[1:-1]:
            raise SystemExit(f"macro defined twice with different values: {name}")
    if body["r@matched"] != body["hl@matched"]:
        raise SystemExit("the ratio table and the campaigns count different matched pairs")
    # Short names keep scale and record keys out of the source of the main text.
    for alias, name in (("hl@agree@chill", "a@liquid@1e3@dseams-chill@ovito-chill@pct"),
                        ("hl@agree@q6", "a@liquid@1e3@dseams-q6@freud-q6@pct"),
                        ("hl@timings", "h@records"),
                        ("hl@cycless@six", "k@confined@1e3@genice-rings@six"),
                        ("hl@dseams@six", "k@confined@1e3@dseams-rings@six")):
        define(alias, body[name])
    header = ("% Generated by stats.py from the records in this directory.\n"
              "% snapper:off\n"
              "\\providecommand{\\dsv}[1]{\\csname dsv@#1\\endcsname}\n")
    (DATA / "macros.tex").write_text(header + "\n".join(macros) + "\n% snapper:on\n")
    sys.stdout.write(f"{len(macros)} macros\n")


if __name__ == "__main__":
    main()
