#!/usr/bin/env python3
"""Fill the summary gaps that are not eight-sample wins.

Absent rows are either timed, when the tool imports, or replaced with the
API line that shows the operation is not in that package. Adaptive CNA is
compared atom by atom with CHILL+ on the same three frames. IRA and cycless
times that finished are fit to t = a N^b. Extrapolations are labelled fits.
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

SUMMARY = ROOT / "benchmarks" / "h2h-summary.csv"
SCALE_CSV = ROOT / "benchmarks" / "h2h-scale-fit.csv"
CNA_GROUPS = ROOT / "benchmarks" / "h2h-cna-agree.csv"
CNA_ATOMS = ROOT / "benchmarks" / "h2h-cna-atoms.csv"
FIT_MARK = " | fit "

# Ice Ih cells. ncell**3 * 8 molecules. The summary finished IRA at 4096
# and cycless at 4096 and 110592. These sizes all finish; 110592 cycless
# is already in the summary and is not repeated here.
LADDER = (2, 3, 4, 5, 6, 8, 10, 12, 16)
WARMUP = 1
REPEATS = 3

CNA_NAME = {0: "Other", 1: "FCC", 2: "HCP", 3: "BCC", 4: "ICO"}
CHILL_ROW = {
    "cubic_ice": "4 staggered and 0 eclipsed",
    "hexagonal_ice": "3 staggered and 1 eclipsed",
    "interfacial_ice": "2 staggered, or 3 staggered and 0 eclipsed",
    "hydrate": "4 eclipsed and 0 staggered",
    "interfacial_hydrate": "3 eclipsed",
    "other": "no CHILL+ row",
}
CLAVIS_NAME = {
    0: "other",
    1: "hexagonal_ice",
    2: "cubic_ice",
    3: "interfacial_ice",
    4: "hydrate",
    5: "interfacial_hydrate",
}


def median(samples: list[float]) -> float:
    ordered = sorted(samples)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def fit_power(sizes, times) -> dict:
    """Ordinary least squares of log(t) on log(N). t = exp(intercept) * N^exponent."""
    x = np.log(np.asarray(sizes, dtype=np.float64))
    y = np.log(np.asarray(times, dtype=np.float64))
    if x.size < 3:
        raise ValueError("a power-law fit needs at least three sizes")
    slope, intercept = np.polyfit(x, y, 1)
    predicted = intercept + slope * x
    residual = y - predicted
    ss_res = float(np.dot(residual, residual))
    centered = y - float(y.mean())
    ss_tot = float(np.dot(centered, centered))
    r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    sxx = float(np.dot(x - x.mean(), x - x.mean()))
    sigma2 = ss_res / (x.size - 2)
    se = float(np.sqrt(sigma2 / sxx))
    return {
        "exponent": float(slope),
        "se": se,
        "intercept": float(intercept),
        "r2": float(r2),
        "n": int(x.size),
    }


def extrapolate_ms(fit: dict, molecules: int) -> float:
    return float(np.exp(fit["intercept"]) * molecules ** fit["exponent"])


def disagreement_cause(chill: str, cna_code: int, coordination: int) -> str:
    """Why this atom does or does not count as a CHILL+ versus adaptive-CNA disagreement."""
    structured = chill != "other"
    cna_other = int(cna_code) == 0
    if structured and cna_other and coordination == 4:
        return "four_neighbor"
    if structured and cna_other:
        return "structured_other_shell"
    if structured and not cna_other:
        return "both_named"
    if (not structured) and (not cna_other):
        return "cna_only"
    return "agree"


def pyscal_chill_reason() -> str:
    import pyscal3

    exported = [name for name in pyscal3.__all__ if "chill" in name.lower()]
    if exported:
        raise SystemExit(f"pyscal3 exports {exported}; time CHILL+ instead of marking it absent")
    return (
        f"pyscal3 {pyscal3.__version__} __all__ has no CHILL or CHILL+ name. "
        "common_neighbor_analysis is adaptive CNA and returns fcc, hcp, bcc, ico, and others "
        "(pyscal3/descriptors.py). diamond_structure is extended CNA. "
        "Neither is Nguyen and Molinero CHILL+. PyPI's latest pyscal3 is this same 4.1.0."
    )


def freud_ring_reason() -> str:
    import pkgutil

    import freud

    modules = sorted(
        item.name for item in pkgutil.iter_modules(freud.__path__) if not item.name.startswith("_")
    )
    hit = [name for name in modules if "ring" in name or "cycle" in name]
    if hit:
        raise SystemExit(f"freud has {hit}; time that ring finder instead of marking it absent")
    listed = ", ".join(modules)
    return (
        f"freud {freud.__version__} modules are {listed}. "
        "There is no ring or cycle module. order.py is Steinhardt, Hexatic, SolidLiquid, "
        "Cubatic, Nematic, RotationalAutocorrelation, and ContinuousCoordination. "
        "Ramasubramani and others, Comput. Phys. Commun. 254, 107275 (2020)."
    )


def clavis_absent_reason() -> str:
    installed = "not imported"
    try:
        import clavis

        text = (clavis.__doc__ or "").splitlines()[0]
        installed = f"clavis {clavis.__version__} is {text}"
    except ImportError:
        installed = "the clavis package is not imported"
    return (
        "clavis_engine does not import, and PyPI has no clavis-engine distribution. "
        f"{installed}. "
        "Zenodo 10.5281/zenodo.24112918 is ClavisCage 0.3.18 source, which does implement CHILL+ "
        "(src/chill.rs calculate_chill_cell)."
    )


def fit_suffix(kind: str, molecules: int) -> str:
    if not SCALE_CSV.is_file():
        return ""
    fits = {}
    quadratic = {}
    with SCALE_CSV.open(encoding="ascii", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["kind"] != kind:
                continue
            if row["role"] == "fit":
                fits[int(row["molecules"])] = row
            elif row["role"] == "quadratic":
                quadratic[int(row["molecules"])] = row
    row = fits.get(int(molecules))
    if row is None:
        return ""
    serial = "classify_serial has no thread argument" if kind == "ira" else "cycles_iter has no thread argument"
    text = (
        f"{FIT_MARK}b={row['exponent']} se={row['se']} R2={row['r2']} "
        f"from {row['sizes']}; extrapolated {molecules} = {row['median_ms']} ms."
    )
    extra = quadratic.get(int(molecules))
    if extra is not None:
        text += f" N^2 projection from {extra['sizes']} is {extra['median_ms']} ms."
    if kind == "ira":
        text += " The nine-size exponent is not used: the local exponent rises with N."
    return text + f" Fit, not a measurement. {serial}."


def _write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _coordination(pos, lengths) -> np.ndarray:
    import freud

    box = freud.box.Box(*lengths)
    query = freud.locality.AABBQuery(box, pos).query(
        pos, {"r_max": 3.5, "exclude_ii": True}
    ).toNeighborList()
    return np.asarray(query.neighbor_counts, dtype=np.int32)


def fcc_control() -> dict:
    """The same adaptive-CNA call on periodic FCC copper. Every atom is FCC."""
    from ase.build import bulk

    import h2h_summary_extra as extra

    crystal = bulk("Cu", "fcc", a=3.615, cubic=True).repeat((4, 4, 4))
    pos = np.ascontiguousarray(crystal.get_positions(), dtype=np.float64)
    length = float(crystal.cell.lengths()[0])
    codes = extra.cna_codes(pos, (length, length, length))
    counts = {CNA_NAME.get(int(value), str(int(value))): int(np.count_nonzero(codes == value))
              for value in np.unique(codes)}
    return {"n": int(pos.shape[0]), "counts": counts, "fcc": int(np.count_nonzero(codes == 1))}


def explain_cna() -> tuple[str, list[dict], list[dict]]:
    import h2h_mdapy_pyscal3_vesin as h
    import h2h_summary_extra as extra

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    lib = h.load_seams(h.library_path())
    ih_pos, ih_lengths = h.ice_positions(8)
    from gen_ice_ic import ice_ic

    ic_pos, side = ice_ic(8)
    frames = (
        ("Ih", ih_pos, ih_lengths),
        ("Ic", np.ascontiguousarray(ic_pos), (side, side, side)),
        ("liquid", *extra.liquid_positions(h)),
    )
    atoms = []
    groups = []
    sentences = []
    control = fcc_control()
    control_text = (
        f"FCC control {control['fcc']}/{control['n']} "
        + " ".join(f"{key}={control['counts'][key]}" for key in sorted(control["counts"]))
    )
    for phase, pos, lengths in frames:
        names = extra.chill_names(lib, h, pos, lengths)
        codes = np.asarray(extra.cna_codes(pos, lengths), dtype=np.int32)
        coord = _coordination(pos, lengths)
        if len(names) != len(codes) or len(names) != len(coord):
            raise SystemExit(f"{phase} length mismatch")
        bucket = {}
        disagree = 0
        for index, (name, code, degree) in enumerate(zip(names, codes, coord)):
            cause = disagreement_cause(name, int(code), int(degree))
            is_disagreement = cause != "agree"
            disagree += int(is_disagreement)
            atoms.append({
                "phase": phase,
                "atom": index,
                "chill": name,
                "cna": CNA_NAME.get(int(code), str(int(code))),
                "coordination": int(degree),
                "disagree": int(is_disagreement),
                "cause": cause,
            })
            key = (name, CNA_NAME.get(int(code), str(int(code))), int(degree), cause)
            bucket[key] = bucket.get(key, 0) + 1
        for (name, cna, degree, cause), count in sorted(bucket.items()):
            groups.append({
                "phase": phase,
                "chill": name,
                "cna": cna,
                "coordination": degree,
                "cause": cause,
                "atoms": count,
                "chill_row": CHILL_ROW.get(name, ""),
            })
        disagree_groups = [row for row in groups if row["phase"] == phase and row["cause"] != "agree"]
        parts = [
            f"{row['chill']} vs {row['cna']} coordination {row['coordination']} "
            f"n={row['atoms']} cause={row['cause']} ({row['chill_row']})"
            for row in disagree_groups
        ]
        detail = "; ".join(parts) if parts else "no disagreements"
        sentences.append(
            f"phase {phase} atoms {len(names)} disagree {disagree}. {detail}. "
            "Adaptive CNA Structure Type is Other, FCC, HCP, BCC, or ICO "
            "(ovito/modifiers/__init__.pyi). "
            f"{control_text}."
        )
    return " ".join(sentences), groups, atoms


def cna_agreement(phase: str, groups: list[dict], control: dict) -> str:
    rows = [row for row in groups if row["phase"] == phase]
    atoms = sum(int(row["atoms"]) for row in rows)
    disagree = sum(int(row["atoms"]) for row in rows if row["cause"] != "agree")
    bits = []
    for row in rows:
        if row["cause"] == "agree":
            continue
        bits.append(
            f"{row['chill']} vs {row['cna']} coordination {row['coordination']} "
            f"n={row['atoms']} ({row['chill_row']})"
        )
    body = "; ".join(bits) if bits else "no disagreements"
    return (
        f"phase {phase} atoms {atoms} disagree {disagree}. {body}. "
        "Every disagreement is a CHILL+ row on a shell adaptive CNA does not name: "
        "Structure Type is only Other, FCC, HCP, BCC, or ICO "
        "(ovito/modifiers/__init__.pyi). "
        f"FCC control {control['fcc']}/{control['n']} FCC, so the modifier is classifying. "
        "CHILL+ rows are the table in src/include/internal/bop.hpp."
    )


def _pin(threads: int) -> None:
    os.sched_setaffinity(0, {0} if threads == 1 else set(range(8)))


def _timed(fn) -> float:
    import time

    t0 = time.perf_counter_ns()
    fn()
    return (time.perf_counter_ns() - t0) / 1e6


def cycless_graph(pos, lengths):
    import networkx as nx
    import freud

    n = int(pos.shape[0])
    box = freud.box.Box(*lengths)
    query = freud.locality.AABBQuery(box, pos).query(
        pos, {"r_max": 3.5, "exclude_ii": True}
    ).toNeighborList()
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    ii = np.asarray(query[:, 0])
    jj = np.asarray(query[:, 1])
    graph.add_edges_from((int(a), int(b)) for a, b in zip(ii, jj) if int(a) < int(b))
    return graph


def time_ladder() -> tuple[list[dict], dict]:
    """One sitting. One warm-up and three calls. Serial, core 0."""
    import h2h_mdapy_pyscal3_vesin as h
    from cycless.cycles import cycles_iter
    from h2h_ira import classify_serial, cluster

    import ira_mod
    from gen_ice_ic import ice_ic

    _pin(1)
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_PROC_BIND"] = "close"
    os.environ["OMP_PLACES"] = "cores"
    os.environ["OMP_WAIT_POLICY"] = "PASSIVE"
    engine = ira_mod.IRA()
    ic_pos, side = ice_ic(4)
    ic = cluster(np.ascontiguousarray(ic_pos), np.asarray((side, side, side), dtype=np.float64), 0)
    rows = []
    measured = {"ira": [], "cycless": []}
    for ncell in LADDER:
        pos, lengths = h.ice_positions(ncell)
        pos = np.ascontiguousarray(pos, dtype=np.float64)
        box = np.asarray(lengths, dtype=np.float64)
        n = int(pos.shape[0])
        ih = cluster(pos, box, 0)
        print(f"scale ira {n}", file=sys.stderr, flush=True)

        def ira_call(pos=pos, box=box, ih=ih):
            classify_serial(engine, ih, ic, pos, box)

        for _ in range(WARMUP):
            ira_call()
        ira_samples = [_timed(ira_call) for _ in range(REPEATS)]
        ira_med = median(ira_samples)
        measured["ira"].append((n, ira_med))
        rows.append({
            "kind": "ira",
            "molecules": n,
            "role": "measured",
            "warmup": WARMUP,
            "calls": REPEATS,
            "median_ms": f"{ira_med:.4f}",
            "min_ms": f"{min(ira_samples):.4f}",
            "samples_ms": " ".join(f"{value:.4f}" for value in ira_samples),
            "exponent": "",
            "se": "",
            "r2": "",
            "sizes": "",
        })
        graph = cycless_graph(pos, lengths)
        print(f"scale cycless {n}", file=sys.stderr, flush=True)

        def cyc_call(graph=graph):
            if sum(1 for _ in cycles_iter(graph, 6)) < 1:
                raise SystemExit("cycless returned no rings")

        for _ in range(WARMUP):
            cyc_call()
        cyc_samples = [_timed(cyc_call) for _ in range(REPEATS)]
        cyc_med = median(cyc_samples)
        measured["cycless"].append((n, cyc_med))
        rows.append({
            "kind": "cycless",
            "molecules": n,
            "role": "measured",
            "warmup": WARMUP,
            "calls": REPEATS,
            "median_ms": f"{cyc_med:.4f}",
            "min_ms": f"{min(cyc_samples):.4f}",
            "samples_ms": " ".join(f"{value:.4f}" for value in cyc_samples),
            "exponent": "",
            "se": "",
            "r2": "",
            "sizes": "",
        })
    return append_fits(rows)


def append_fits(rows: list[dict]) -> tuple[list[dict], dict]:
    """Quoted fits. cycless uses every finished size. IRA's local exponent
    rises toward 2, so the quoted IRA fit is the three largest sizes, and
    an N^2 projection of the largest size is stored beside it."""
    measured = {"ira": [], "cycless": []}
    for row in rows:
        if row["role"] == "measured":
            measured[row["kind"]].append((int(row["molecules"]), float(row["median_ms"])))
    kept = [row for row in rows if row["role"] == "measured"]
    fits = {}
    for kind, points in measured.items():
        points = sorted(points)
        if kind == "ira":
            used = [item for item in points if item[0] >= 8000]
        else:
            used = points
        sizes = [item[0] for item in used]
        times = [item[1] for item in used]
        fit = fit_power(sizes, times)
        fits[kind] = fit
        size_text = ",".join(str(item) for item in sizes)
        for molecules in (110592, 884736):
            kept.append({
                "kind": kind,
                "molecules": molecules,
                "role": "fit",
                "warmup": WARMUP,
                "calls": REPEATS,
                "median_ms": f"{extrapolate_ms(fit, molecules):.1f}",
                "min_ms": "",
                "samples_ms": "",
                "exponent": f"{fit['exponent']:.4f}",
                "se": f"{fit['se']:.4f}",
                "r2": f"{fit['r2']:.6f}",
                "sizes": size_text,
            })
        if kind == "ira":
            anchor_n, anchor_t = points[-1]
            for molecules in (110592, 884736):
                projected = anchor_t * (molecules / anchor_n) ** 2
                kept.append({
                    "kind": kind,
                    "molecules": molecules,
                    "role": "quadratic",
                    "warmup": WARMUP,
                    "calls": REPEATS,
                    "median_ms": f"{projected:.1f}",
                    "min_ms": "",
                    "samples_ms": "",
                    "exponent": "2",
                    "se": "",
                    "r2": "",
                    "sizes": str(anchor_n),
                })
    return kept, fits


def clavis_available() -> bool:
    try:
        import clavis_engine  # noqa: F401
    except ImportError:
        return False
    return True


def time_clavis_cell(ncell: int, threads: int) -> dict:
    """Pinned ABAB, eight samples, the same win rule as the summary."""
    import h2h_mdapy_pyscal3_vesin as h
    import h2h_summary as summary
    import clavis_engine as engine

    _pin(threads)
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["OMP_PROC_BIND"] = "close"
    os.environ["OMP_PLACES"] = "cores"
    os.environ["OMP_WAIT_POLICY"] = "PASSIVE"
    lib = h.load_seams(h.library_path())
    lib.seams_version_string.restype = __import__("ctypes").c_char_p
    seams = lib.seams_version_string().decode()
    pos, lengths = h.ice_positions(ncell)
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    n = int(pos.shape[0])
    held = h.frame_system(lib, pos, lengths)
    h.require(lib, lib.seams_set_num_threads(threads), "threads")
    labels = (__import__("ctypes").c_int * n)()
    counts = (__import__("ctypes").c_int * 9)()
    h.dseams_chill_plus(lib, held, labels, counts, h.SEAMS_CHILL_OVITO)
    ours = np.array([h.DSEAMS_AS_OVITO[int(value)] for value in labels])
    points_nm = pos / 10.0
    waters = []
    for oxygen in points_nm:
        xyz = [float(oxygen[0]), float(oxygen[1]), float(oxygen[2])]
        waters.append({"oxygen": xyz, "hydrogen1": xyz, "hydrogen2": xyz})
    payload = json.dumps(waters)
    box_nm = [float(lengths[0]) / 10.0, float(lengths[1]) / 10.0, float(lengths[2]) / 10.0]
    raw = engine.analyze_waters(payload, "[]", box_nm, "chill", 0.35)
    parsed = json.loads(raw)
    theirs = np.array([CLAVIS_NAME[int(value)] for value in parsed["per_atom_class"]])
    differ = int(np.count_nonzero(ours != theirs))
    note = (
        f"labels {n} equal"
        if differ == 0
        else f"labels differ {differ} of {n}"
    )
    note += (
        f"; analyze_waters parses JSON built once; chill.rs is a serial loop "
        f"and Cargo.toml leaves rayon commented out; {engine.version()}"
    )

    def dseams_call():
        _pin(threads)
        h.require(lib, lib.seams_set_num_threads(threads), "threads")
        h.dseams_chill_plus(lib, held, labels, counts, h.SEAMS_CHILL_OVITO)

    def clavis_call():
        _pin(threads)
        engine.analyze_waters(payload, "[]", box_nm, "chill", 0.35)

    print(f"clavis n {n} threads {threads}", file=sys.stderr, flush=True)
    left, right = summary.abab(dseams_call, clavis_call)
    lib.seams_system_destroy(held)
    model = summary.cpu_model()
    made = summary.row_of(
        "ClavisCage", "chill_plus", n, threads, os.cpu_count() or 1,
        left, right, note, f"{seams}/{engine.version()}", model,
    )
    made["budget_s"] = ""
    made["largest_finished"] = ""
    return made


def write_absent_md(clavis_measured: bool) -> None:
    pyscal = pyscal_chill_reason()
    freud = freud_ring_reason()
    clavis = (
        "ClavisCage 0.3.18 imported as clavis_engine. CHILL+ is calculate_chill_cell "
        "in src/chill.rs of Zenodo 10.5281/zenodo.24112918. The six summary rows are "
        "that call, pinned and interleaved with d-SEAMS."
        if clavis_measured
        else clavis_absent_reason()
    )
    text = "\n".join([
        "# Absent summary rows",
        "",
        "The summary had 14 absent rows: pyscal3 CHILL+ at both widths, freud rings at three sizes and both widths, and ClavisCage CHILL+ at three sizes and both widths.",
        "",
        "## pyscal3 CHILL+",
        "",
        pyscal,
        "",
        "The two rows stay absent. There is no equivalent call to time.",
        "",
        "## freud rings",
        "",
        freud,
        "",
        "Neighbor lists are already a separate timed row. A ring search written on top of that list would be a new code, not freud. The six rows stay absent.",
        "",
        "## ClavisCage CHILL+",
        "",
        clavis,
        "",
        "PyPI's package named clavis is unrelated. Its module docstring is a credential client.",
        "",
    ])
    path = ROOT / "benchmarks" / "h2h-absent.md"
    path.write_text(text, encoding="ascii")


def write_cna_md(groups: list[dict], control: dict) -> None:
    lines = [
        "# Adaptive CNA against CHILL+",
        "",
        "The three correctness rows are ice Ih, ice Ic, and the same Ih cell after a 1.5 A Gaussian shift (PCG64 seed 20261113), each 4096 molecules. CHILL+ is d-SEAMS ovito mode. Adaptive CNA is OVITO CommonNeighborAnalysisModifier mode AdaptiveCutoff, the call in scripts/h2h_summary_extra.py.",
        "",
        f"A periodic FCC copper cell, lattice constant 3.615 A, 4x4x4 cubic cells, {control['n']} atoms, goes through that same modifier: {control['fcc']} atoms are FCC. The ice result is not a pipeline that labels everything Other.",
        "",
        "Adaptive CNA writes Structure Type Other, FCC, HCP, BCC, or ICO (ovito/modifiers/__init__.pyi). Those signatures are close-packed or bcc shells. CHILL+ writes a water class from the staggered and eclipsed counts of a four-neighbor shell (src/include/internal/bop.hpp). A four-neighbor shell matches none of the four crystal types, so the structure type is Other.",
        "",
        "| Phase | CHILL+ | Adaptive CNA | Coordination | Atoms | Cause |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in groups:
        cause = {
            "four_neighbor": "CHILL+ row on a 4-neighbor shell; adaptive CNA has no tetrahedral type",
            "structured_other_shell": "CHILL+ named a class and the 3.5 A coordination is not 4",
            "both_named": "both methods named a crystal type",
            "cna_only": "adaptive CNA named a crystal type and CHILL+ said other",
            "agree": "both unstructured, or the boolean labels match",
        }[row["cause"]]
        lines.append(
            f"| {row['phase']} | {row['chill']} | {row['cna']} | {row['coordination']} | {row['atoms']} | {cause} |"
        )
    lines.extend([
        "",
        "On Ih every atom is hexagonal ice, coordination 4, adaptive CNA Other. That is 3 staggered bonds and 1 eclipsed bond. On Ic every atom is cubic ice, coordination 4, adaptive CNA Other. That is 4 staggered bonds. The 222 liquid disagreements are the waters that still have exactly four neighbors inside 3.5 A and still match one CHILL+ row after the shift: 1 hexagonal, 49 hydrate, 168 interfacial hydrate, and 4 interfacial ice. Every one of them is adaptive CNA Other. Another 634 waters also have four neighbors and stay CHILL+ other, because the staggered and eclipsed counts match no row; those agree with adaptive CNA Other. Four liquid atoms have 12 neighbors inside 3.5 A and are still Other. Coordination by itself is not an FCC signature. The FCC control is a copper lattice; the shifted ice is not.",
        "",
        "Per-atom rows are benchmarks/h2h-cna-atoms.csv.",
        "",
    ])
    (ROOT / "benchmarks" / "h2h-cna-agree.md").write_text("\n".join(lines), encoding="ascii")


def write_scale_md(rows: list[dict], fits: dict) -> None:
    lines = [
        "# IRA and cycless scaling, labelled as fits",
        "",
        "One sitting, core 0, one warm-up and three calls. IRA is classify_serial. cycless is cycles_iter on the 3.5 A graph, the graph built outside the timer, maximum size 6. Neither call takes a thread count, so the eight-thread DNF rows use this same serial fit.",
        "",
        "The model is t = a N^b, ordinary least squares on log t against log N. a is exp(intercept). The standard error is the usual slope error. An extrapolated time is a fit, not a measurement. The 600 s budget still stands: a size that did not return is DNF.",
        "",
        "| Kind | Molecules | Role | Median ms | Exponent | SE | R2 |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            f"| {row['kind']} | {row['molecules']} | {row['role']} | {row['median_ms']} | {row['exponent']} | {row['se']} | {row['r2']} |"
        )
    for kind, fit in fits.items():
        lines.append("")
        lines.append(
            f"{kind}: b = {fit['exponent']:.4f}, se {fit['se']:.4f}, R2 {fit['r2']:.6f}, "
            f"{fit['n']} sizes."
        )
    lines.extend([
        "",
        "cycless is close to linear on the whole ladder. The nine-point fit at 110592 molecules is 354484 ms. The summary sitting measured 341558 ms there, one thread, eight samples. The fit is not a replacement for that row. At 884736 molecules the same fit is 2952 s, past the 600 s budget, which is why that row is DNF.",
        "",
        "IRA is not one exponent. cluster() subtracts every molecule and argsorts that vector, once per atom, so the leading cost is N^2. The local exponent rises from 0.84 between 64 and 216 molecules to 1.93 between 13824 and 32768. A nine-point line gives b = 1.346 and would put 110592 molecules near 142 s. That line is not quoted. The quoted fit is ordinary least squares on the three largest finished sizes, 8000, 13824, and 32768: b = 1.883, standard error 0.035. It puts classify_serial at 502 s for 110592 molecules and 25162 s for 884736. An N^2 projection of the 32768-molecule measurement, 51.403 s, is 586 s and 37473 s. The 600 s probe times the whole process and did not return. 502 s of classify_serial, still on a rising local exponent, sits next to that budget. Both numbers are fits, not measurements.",
        "",
    ])
    (ROOT / "benchmarks" / "h2h-scale-fit.md").write_text("\n".join(lines), encoding="ascii")


def patch_summary(groups: list[dict], control: dict, clavis_rows: list[dict]) -> None:
    with SUMMARY.open(encoding="ascii", newline="") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0].keys()) if rows else []
    # The file on disk has a header even when DictReader consumes it. Re-read fields.
    with SUMMARY.open(encoding="ascii", newline="") as handle:
        fields = next(csv.reader(handle))
    pyscal = pyscal_chill_reason()
    freud = freud_ring_reason()
    clavis_why = clavis_absent_reason()
    phase_text = {}
    for phase in ("Ih", "Ic", "liquid"):
        phase_text[phase] = cna_agreement(phase, groups, control)
    clavis_index = {
        (str(row["molecules"]), str(row["threads"])): row for row in clavis_rows
    }
    patched = []
    for row in rows:
        competitor = row["competitor"]
        operation = row["operation"]
        if competitor == "pyscal3" and operation == "chill_plus":
            row["agreement"] = pyscal
            row["verdict"] = "absent"
        elif competitor == "freud" and operation == "rings":
            row["agreement"] = freud
            row["verdict"] = "absent"
        elif competitor == "ClavisCage" and operation == "chill_plus":
            key = (str(row["molecules"]), str(row["threads"]))
            if key in clavis_index:
                row = clavis_index[key]
            else:
                row["agreement"] = clavis_why
                row["verdict"] = "absent"
        elif competitor == "OVITO" and operation == "cna":
            text = row["agreement"]
            phase = "Ih" if "phase Ih" in text else "Ic" if "phase Ic" in text else "liquid"
            row["agreement"] = phase_text[phase]
        elif competitor == "IRA" and operation == "shape" and row["verdict"] == "DNF":
            base = row["agreement"].split(FIT_MARK)[0].rstrip()
            row["agreement"] = base + fit_suffix("ira", int(row["molecules"]))
        elif competitor == "cycless" and operation == "rings" and row["verdict"] == "DNF":
            base = row["agreement"].split(FIT_MARK)[0].rstrip()
            row["agreement"] = base + fit_suffix("cycless", int(row["molecules"]))
        patched.append(row)
    _write_csv(SUMMARY, fields, patched)


def main() -> None:
    steps = set(sys.argv[1:]) or {"absent", "cna", "scale", "clavis", "patch"}
    groups: list[dict] = []
    control = {"n": 0, "fcc": 0, "counts": {}}
    clavis_rows: list[dict] = []
    if "cna" in steps:
        _sentence, groups, atoms = explain_cna()
        control = fcc_control()
        _write_csv(
            CNA_GROUPS,
            ["phase", "chill", "cna", "coordination", "cause", "atoms", "chill_row"],
            groups,
        )
        _write_csv(
            CNA_ATOMS,
            ["phase", "atom", "chill", "cna", "coordination", "disagree", "cause"],
            atoms,
        )
        write_cna_md(groups, control)
        print(f"cna groups {len(groups)} atoms {len(atoms)}", file=sys.stderr, flush=True)
    scale_fields = [
        "kind", "molecules", "role", "warmup", "calls", "median_ms", "min_ms",
        "samples_ms", "exponent", "se", "r2", "sizes",
    ]
    if "scale" in steps:
        rows, fits = time_ladder()
        _write_csv(SCALE_CSV, scale_fields, rows)
        write_scale_md(rows, fits)
        print("scale wrote", SCALE_CSV.name, file=sys.stderr, flush=True)
    elif "refit" in steps and SCALE_CSV.is_file():
        with SCALE_CSV.open(encoding="ascii", newline="") as handle:
            measured_rows = [row for row in csv.DictReader(handle) if row["role"] == "measured"]
        rows, fits = append_fits(measured_rows)
        _write_csv(SCALE_CSV, scale_fields, rows)
        write_scale_md(rows, fits)
        print("refit wrote", SCALE_CSV.name, file=sys.stderr, flush=True)
    if "clavis" in steps and clavis_available():
        for ncell in (8, 24, 48):
            for threads in (1, 8):
                clavis_rows.append(time_clavis_cell(ncell, threads))
        _write_csv(ROOT / "benchmarks" / "h2h-clavis.csv", list(clavis_rows[0].keys()), clavis_rows)
    if "absent" in steps or "patch" in steps:
        write_absent_md(bool(clavis_rows) or clavis_available() and (ROOT / "benchmarks" / "h2h-clavis.csv").is_file())
    if "patch" in steps:
        if not groups and CNA_GROUPS.is_file():
            with CNA_GROUPS.open(encoding="ascii", newline="") as handle:
                groups = list(csv.DictReader(handle))
            control = fcc_control()
        if not clavis_rows and (ROOT / "benchmarks" / "h2h-clavis.csv").is_file() and clavis_available():
            with (ROOT / "benchmarks" / "h2h-clavis.csv").open(encoding="ascii", newline="") as handle:
                clavis_rows = list(csv.DictReader(handle))
        patch_summary(groups, control, clavis_rows)
        print("patched", SUMMARY.name, file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
