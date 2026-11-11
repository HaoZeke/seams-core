# Copyright (c) d-SEAMS developers
"""Regenerate summary.tsv and CLAIMS.tsv from the raw records.

Median and quartiles are Hyndman-Fan type 7. The interquartile range is
the type 7 third quartile minus the type 7 first quartile.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RECORDS, parse_counts, parse_repeats, summary

ROOT_RECORDS = RECORDS


def load_text(path: Path) -> str:
    return path.read_text() if path.is_file() else ""


def series_rows() -> list[str]:
    rows = ["kind\tmeasurement\tsystem\tquantity\tn\tmedian\tq1\tq3\tiqr\tunit\tsource"]
    for path in sorted(ROOT_RECORDS.rglob("*")):
        if path.suffix not in {".csv", ".log", ".txt"}:
            continue
        text = path.read_text(errors="replace")
        if "repeat " not in text:
            continue
        repeats = parse_repeats(text)
        rel = path.relative_to(ROOT_RECORDS).as_posix()
        measurement = rel.split("/", 1)[0]
        system = path.stem
        for quantity, values in sorted(repeats.items()):
            if quantity in {"largest", "ice_like", "subset"}:
                stat = summary(values)
                rows.append(
                    "count\t%s\t%s\t%s\t%d\t%.6g\t%.6g\t%.6g\t%.6g\tcount\t%s"
                    % (measurement, system, quantity, len(values), stat["median"],
                       stat["q1"], stat["q3"], stat["iqr"], rel)
                )
                continue
            if not quantity.endswith("_ms"):
                continue
            stat = summary(values)
            rows.append(
                "time\t%s\t%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\tms\t%s"
                % (measurement, system, quantity, len(values), stat["median"],
                   stat["q1"], stat["q3"], stat["iqr"], rel)
            )
        if path.suffix == ".csv" and "count " in text:
            for key, value in parse_counts(text).items():
                rows.append(
                    "count\t%s\t%s\t%s\t1\t%s\t\t\t\tcount\t%s"
                    % (measurement, system, key, value, rel)
                )
    return rows


def cell(path: Path, quantity: str) -> str:
    text = load_text(path)
    values = parse_repeats(text).get(quantity, [])
    if not values:
        return ""
    stat = summary(values)
    return "%.6f (q1 %.6f, q3 %.6f, iqr %.6f, n %d)" % (
        stat["median"], stat["q1"], stat["q3"], stat["iqr"], len(values)
    )


def overall(path: Path) -> str:
    if not path.is_file():
        return ""
    for line in path.read_text().splitlines():
        if line.startswith("overall"):
            parts = line.split()
            return parts[4] if len(parts) > 4 else ""
    return ""


def count_of(path: Path, key: str) -> str:
    return parse_counts(load_text(path)).get(key, "")


def claims() -> list[str]:
    b = ROOT_RECORDS / "b-cages"
    a = ROOT_RECORDS / "a-chill"
    rows = [
        "paper_sentence\tpaper_value\tnew_value\tsource_line\tcommand\tstatus",
        (
            "On a 1000-molecule ice Ic frame the whole path takes 5.28 ms.\t"
            "5.28 ms\t%s\t"
            "paper/tex/rg_main.tex:508\t"
            "measure_cpc time cages benchmarks/cpc/frames/ic-1000.lammpstrj 1 1 1 6\t"
            "NOT_RERUN"
        ) % cell(b / "ic-1000.csv", "path_ms"),
        (
            "The hexagonal-cage stage on that frame takes 1.24 ms.\t"
            "1.24 ms\t%s\t"
            "paper/tex/figures/data/macros.tex:306\t"
            "measure_cpc time cages benchmarks/cpc/frames/ic-1000.lammpstrj 1 1 1 6\t"
            "NOT_RERUN"
        ) % cell(b / "ic-1000.csv", "hc_ms"),
        (
            "The double-diamond stage, packing included, takes 2.61 ms.\t"
            "2.61 ms\t%s\t"
            "paper/tex/figures/data/macros.tex:314\t"
            "measure_cpc time cages benchmarks/cpc/frames/ic-1000.lammpstrj 1 1 1 6\t"
            "slowdown; cage counts identical"
        ) % cell(b / "ic-1000.csv", "ddc_ms"),
        (
            "Packing turns 2000 equators into 500 double-diamond cages on 1000 molecules of ice Ic.\t"
            "2000 equators, 500 cages\t"
            "equators %s cages %s\t"
            "paper/tex/rg_supplement.tex:278\t"
            "measure_cpc time cages benchmarks/cpc/frames/ic-1000.lammpstrj 1 1 1 6\t"
            "NOT_RERUN"
        ) % (count_of(b / "ic-1000.csv", "ddc_equators"), count_of(b / "ic-1000.csv", "ddc_cages")),
        (
            "Every molecule of ideal ice Ih lies on hexagonal cages: 500 HC and 0 DDC at 1000 molecules.\t"
            "500 HC, 0 DDC\t"
            "hc %s ddc %s\t"
            "paper/tex/figures/data/macros.tex:4301\t"
            "measure_cpc time cages benchmarks/cpc/frames/ih-1000.lammpstrj 1 1 1 6\t"
            "NOT_RERUN"
        ) % (count_of(b / "ih-1000.csv", "hc_cages"), count_of(b / "ih-1000.csv", "ddc_cages")),
        (
            "Fitted exponents of the hexagonal and double-diamond stages on ideal ice Ic are 0.889 and 0.76.\t"
            "0.889 and 0.76\t"
            "current code ddc exponent 0.823 on 216 to 4096 including 2744; 0.826 without 2744\t"
            "paper/tex/figures/data/macros.tex:706\t"
            "measure_cpc time cages benchmarks/cpc/frames/ic-*.lammpstrj\t"
            "slowdown; cage counts identical"
        ),
        (
            "Structure I holds 2 5^12 cages and 6 5^12 6^2 cages.\t"
            "2 and 6\t"
            "see genice_sI-signature.txt\t"
            "paper/tex/figures/data/macros.tex:4310\t"
            "seams cages --graph cutoff --signature 512+51262+51264 input/traj/genice_sI.lammpstrj\t"
            "NOT_RERUN"
        ),
        (
            "On ice Ic and ice Ih, d-SEAMS CHILL+ and OVITO CHILL+ give the same class to every molecule.\t"
            "100 percent\t"
            "ice Ic %s ice Ih %s\t"
            "paper/tex/rg_supplement.tex:650\t"
            "measure_cpc labels chillplus plus ovito ChillPlusModifier cutoff 3.5\t"
            "NOT_RERUN"
        ) % (overall(a / "ice-Ic-chillplus-vs-ovito-merged.tsv"),
             overall(a / "ice-Ih-chillplus-vs-ovito-merged.tsv")),
        (
            "Published CHILL+ and the ovito-chill+ option against OVITO on ideal ice Ic.\t"
            "option absent on main; present at e8610637\t"
            "published %s ovito-chill+ %s\t"
            "benchmarks/cpc/records/a-chill/ice-Ic-tip-chillplus-vs-ovito-merged.tsv:2\t"
            "seams_chill_plus_rule rule chill+ and ovito-chill+ at e8610637\t"
            "NOT_RERUN"
        ) % (overall(a / "ice-Ic-tip-chillplus-vs-ovito-merged.tsv"),
             overall(a / "ice-Ic-tip-ovito-chillplus-vs-ovito-merged.tsv")),
        (
            "CHILL+ agreement on the 512-molecule mW liquid at 215 K.\t"
            "campaign liquid dump was not rerun\t"
            "main chill+ %s published %s ovito-chill+ %s\t"
            "paper/tex/rg_main.tex:473\t"
            "measure_cpc labels chillplus benchmarks/cpc/frames/mw-liquid-512.lammpstrj\t"
            "PRE_FIX_PENDING_CHILLPLUS_PR2"
        ) % (overall(a / "mW-liquid-T215-chillplus-vs-ovito-merged.tsv"),
             overall(a / "mW-liquid-T215-tip-chillplus-vs-ovito-merged.tsv"),
             overall(a / "mW-liquid-T215-tip-ovito-chillplus-vs-ovito-merged.tsv")),
        (
            "Structure I CHILL+ labels against OVITO, and plain CHILL against the same modifier.\t"
            "not a numbered cell in the agreement table\t"
            "chill+ %s chill %s\t"
            "paper/tex/rg_supplement.tex:666\t"
            "measure_cpc labels chillplus input/traj/genice_sI.lammpstrj\t"
            "NOT_RERUN"
        ) % (overall(a / "sI-chillplus-vs-ovito-merged.tsv"),
             overall(a / "sI-chill-vs-ovito-merged.tsv")),
        (
            "Largest ice-like cluster on the shipped mW cubic trajectory, q6 greater than 0.5 and CHILL.\t"
            "not given as a single integer in the results section\t"
            "4096 of 4096 on every frame, both criteria\t"
            "src/cluster.cpp:386\t"
            "measure_cpc time cluster-q6 input/traj/mW_cubic.lammpstrj and cluster-chill\t"
            "NOT_RERUN"
        ),
        (
            "Largest ice-like cluster on the shipped nucleation first frame.\t"
            "not given for this single frame\t"
            "q6 116 chill 155 of 4096\t"
            "src/cluster.cpp:386\t"
            "measure_cpc time cluster-q6 input/traj/nucleation_f1.lammpstrj\t"
            "NOT_RERUN"
        ),
        (
            "Shipped mW cubic frame, 4096 molecules, CHILL+ against OVITO.\t"
            "4096 of 4096\t"
            "%s\t"
            "docs/orgmode/paper/benchmarks.org:218\t"
            "measure_cpc labels chillplus input/traj/mW_cubic.lammpstrj\t"
            "NOT_RERUN"
        ) % overall(a / "mW-cubic-chillplus-vs-ovito-merged.tsv"),
        (
            "Rodger F4 on exampleTraj, oxygen type 2 and hydrogen type 1.\t"
            "finite when hydrogens are kept\t"
            "mean 0.221753 n 250, identical on all 6 repeats\t"
            "tests/test_order_parameter.cpp:454\t"
            "seams f4 input/traj/exampleTraj.lammpstrj --type 2 --htype 1\t"
            "NOT_RERUN"
        ),
        (
            "Ideal ice Ic at sigma 0 is entirely cubic for the cutoff cage labels and for CHILL+.\t"
            "acc 1.000\t"
            "see e-sota/sota.tsv\t"
            "docs/orgmode/paper/benchmarks.org is silent; scripts/sota_compare.py score\t"
            "seams cages --graph cutoff and measure_cpc labels chillplus\t"
            "NOT_RERUN"
        ),
        (
            "Substitutional ions on analytic ice, seams ions --complete.\t"
            "ion_atlas.py seeded_ice 1.000 at low doping\t"
            "see e-ions/ions.tsv\t"
            "dseams2_repro/scripts/ion_atlas.py\t"
            "seams ions --ion-types 3,4 --type 1 --complete\t"
            "NOT_RERUN"
        ),
        (
            "Ice/brine direct-coexistence columns.\t"
            "stored BRINE summary\t"
            "see e-brine/status.txt\t"
            "dseams2_repro/results/brine/conditions.txt\t"
            "python brine_summary.py results/brine\t"
            "NOT_RERUN"
        ),
        (
            "Host-device ladder from 1000 to 32000 molecules.\t"
            "GPU batch record\t"
            "GPU_REQUIRED\t"
            "docs/orgmode/paper/benchmarks.org:292\t"
            "not run\t"
            "NOT_RERUN"
        ),
        (
            "d-SEAMS 1 against d-SEAMS 2 speedups on the 1000-molecule frame.\t"
            "HC about 34300, whole path about 8120\t"
            "not measured; commit 6e51956f is not in this checkout\t"
            "paper/tex/figures/data/tab-cost.tex:12\t"
            "not run\t"
            "NOT_RERUN"
        ),
        (
            "Ring update with three changed bonds is 59.9 times a full search.\t"
            "59.9\t"
            "not measured; the 148-frame nucleation trajectory is not shipped\t"
            "paper/tex/rg_main.tex:637\t"
            "not run\t"
            "NOT_RERUN"
        ),
        (
            "57 v1 decimals with no raw field in REMOVED.tsv.\t"
            "57 decimals\t"
            "not listed; branch submission/cpc-v3 is not on a reachable remote\t"
            "paper repo branch submission/cpc-v3 REMOVED.tsv\t"
            "not run\t"
            "NOT_RERUN"
        ),
        (
            "Main CHILL+ agrees with OVITO on 0.857422 of the 512-molecule mW liquid.\t"
            "0.857422\t"
            "0.857422 is pre-fix, counted before the shell-count bug fix, and pending the CHILL+ pr2 tip; 40 of 73 miss the shell count and 33 are the q3 definition\t"
            "src/bop.cpp:1118\t"
            "benchmarks/cpc/explain_liquid.py\t"
            "PRE_FIX_PENDING_CHILLPLUS_PR2"
        ),
        (
            "Published CHILL+ agrees with OVITO on 0.936 of that liquid; ovito-chill+ agrees on 1.000.\t"
            "0.936 and 1.000\t"
            "479/512 is the published q3 rule against OVITO full-shell average; ovito-chill+ matches OVITO\t"
            "e8610637:src/include/internal/chill_plus_scheme.hpp:28-36\t"
            "explain_tip chill+ and ovito-chill+\t"
            "q3 definition"
        ),
        (
            "Fixed-code CHILL+ liquid agreement after the shell gate.\t"
            "pending\t"
            "pending the CHILL+ pr2 tip; value not measured here\t"
            "src/bop.cpp:1118\t"
            "not run\t"
            "pending CHILL+ pr2 tip"
        ),
        (
            "Tag v2.12.0 double-diamond stage on 1000 ice Ic molecules.\t"
            "2.61 ms\t"
            "2.590806 ms median, six repeats, one warm-up, this machine\t"
            "paper/tex/figures/data/tab-cost.tex:7\t"
            "v2.12.0 measure time cages ic-1000 1 1 1 6\t"
            "rerun"
        ),
        (
            "Current sources double-diamond stage on the same frame.\t"
            "2.61 ms at v2.12.0\t"
            "3.270127 ms; v2.12.0 rerun 2.590806 ms; 2000 equators and 500 cages on both\t"
            "df2a2be3\t"
            "measure_cpc time cages ic-1000\t"
            "slowdown; cage counts identical"
        ),
        (
            "Double-diamond exponent on the current code.\t"
            "0.76\t"
            "0.823 on 216, 512, 1000, 1728, 2744 and 4096; about 0.82\t"
            "df2a2be3\t"
            "log-log least squares on the six paper sizes\t"
            "slowdown; cage counts identical"
        ),
        (
            "Double-diamond timing on perf/ddc-regression-fix.\t"
            "3.270127 ms and exponent 0.823 on current code\t"
            "about 2.74 ms at 1000 molecules, exponent 0.799, commit 2c0db67f\t"
            "perf/ddc-regression-fix 2c0db67f\t"
            "not run in this tree\t"
            "pending verify"
        ),
    ]
    return rows


def main() -> None:
    summary_path = ROOT_RECORDS / "summary.tsv"
    claims_path = ROOT_RECORDS / "CLAIMS.tsv"
    summary_path.write_text("\n".join(series_rows()) + "\n")
    claims_path.write_text("\n".join(claims()) + "\n")
    print(summary_path.name, claims_path.name)


if __name__ == "__main__":
    main()
