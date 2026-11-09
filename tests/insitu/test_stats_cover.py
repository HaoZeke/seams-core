# Copyright (c) d-SEAMS developers
"""Statement and branch coverage for the paper statistics driver."""

import csv
import decimal
import json
import runpy
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "paper" / "tex" / "figures"))
import stats  # noqa: E402

DATA_SRC = REPO / "paper" / "tex" / "figures" / "data"
REL_SOURCES = (
    "scripts/h2h.py",
    "paper/tex/figures/worked_example.py",
    "paper/tex/figures/nucleation_churn.py",
    "paper/tex/repro/updates.sh",
    "tests/bench_churn_sweep.cpp",
    "scripts/baseline_dseams1.cpp",
    "src/topo_one_dim.cpp",
    "paper/tex/figures/stats.py",
    "paper/tex/repro/build.sh",
    "paper/tex/repro/campaign.sh",
    "paper/tex/repro/versions.sh",
)

stats.macros.clear()
stats.NUMBER_LOG.clear()
stats.plot_calls.clear()
stats.SET.clear()
stats.settings()


def _snap():
    return {
        "DATA": stats.DATA,
        "H2H": stats.H2H,
        "ROOT": stats.ROOT,
        "SOURCE": stats.SOURCE,
        "DRAWS": stats.DRAWS,
        "MIN_RING": stats.MIN_RING,
        "PRISM_SIZES": stats.PRISM_SIZES,
        "THRESHOLD_COPIES": stats.THRESHOLD_COPIES,
        "define": stats.define,
        "SET": dict(stats.SET),
        "macros": list(stats.macros),
        "NUMBER_LOG": list(stats.NUMBER_LOG),
        "plot_calls": {key: list(value) for key, value in stats.plot_calls.items()},
    }


def _restore(saved):
    for key, value in saved.items():
        if key == "SET":
            stats.SET.clear()
            stats.SET.update(value)
        elif key == "macros":
            stats.macros.clear()
            stats.macros.extend(value)
        elif key == "NUMBER_LOG":
            stats.NUMBER_LOG.clear()
            stats.NUMBER_LOG.extend(value)
        elif key == "plot_calls":
            stats.plot_calls.clear()
            stats.plot_calls.update(value)
        else:
            setattr(stats, key, value)


@pytest.fixture
def restore():
    saved = _snap()
    try:
        yield
    finally:
        _restore(saved)


def _plant(dest):
    for rel in REL_SOURCES:
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / rel, target)


def _calls(path, rows, atoms):
    lines = [
        "mode cages",
        "atoms %d" % atoms,
        "warmup 1",
        "reps %d" % rows,
        "columns ring_ms hc_ms ddc_ms path_ms",
    ]
    lines.extend("1.0 2.0 3.0 4.0" for _ in range(rows))
    path.write_text("\n".join(lines) + "\n")


def _err(path, six):
    path.write_text(
        "rings %d six_rings %d hc_rings 0 ddc_rings 0 hc_cages 0 ddc_cages 0\n" % (six, six)
    )


def _version_tree(data, mw_rows=30, mw_atoms=8, mw_six=16, bad_six=None):
    for n, calls in stats.SET["v1"].items():
        _calls(data / ("dseams1-ic-%d.txt" % n), calls[1], n)
        six = 1 if bad_six == ("one", n) else 2 * n
        _err(data / ("dseams1-ic-%d.err" % n), six)
    for n in stats.SET["ladder"]:
        _calls(data / ("dseams2-ic-%d.txt" % n), int(stats.SET["v@repeats"]), n)
        six = 1 if bad_six == ("two", n) else 2 * n
        _err(data / ("dseams2-ic-%d.err" % n), six)
    _calls(data / "dseams2-mw.txt", mw_rows, mw_atoms)
    _err(data / "dseams2-mw.err", mw_six)


def _macro_text():
    return "\n".join(stats.macros)


def test_settings_match_the_scripts_and_the_constants():
    assert stats.SET["cutoff"] == decimal.Decimal("3.5")
    assert stats.SET["bond"] == decimal.Decimal("2.75")
    assert stats.SET["q6"] == decimal.Decimal("0.558")
    assert stats.SET["ptm"] == decimal.Decimal("0.1")
    assert stats.SET["repeats"] == decimal.Decimal("6")
    assert stats.SET["draws"] == decimal.Decimal("10000")
    assert stats.SET["o"]["1e6"] == (1, 10)
    assert stats.SET["o"]["1e3"] == (5, 30)
    assert stats.SET["ladder"][-1] == 4096


def test_commit_names_follow_the_record(restore):
    assert stats.commit_for("paper/tex/figures/data/ring-head2head.json") == stats.COMMITS["ring"]
    assert stats.commit_for("records/ring-h2h/run") == stats.COMMITS["ring"]
    assert stats.commit_for("figures/o@1e3") == stats.COMMITS["ring"]
    assert stats.commit_for("paper/tex/figures/data/neigh-head2head.json") == stats.COMMITS["neigh"]
    assert stats.commit_for("dseams1-ic-216.txt") == stats.COMMITS["dseams1"]
    assert stats.commit_for("paper/tex/figures/stats.py") == stats.COMMITS["measure"]


def test_an_empty_sample_has_no_order_statistic(restore):
    with pytest.raises(IndexError):
        stats.median_interval([])


def test_type7_edges(restore):
    with pytest.raises(SystemExit, match="n >= 6"):
        stats.type7_quantile([1, 2, 3, 4, 5], 0.5)
    assert stats.type7_quantile([1, 2, 3, 4, 5, 9], 1) == 9
    assert stats.type7_quantile([1, 2, 3, 4, 5, 9], -0.2) == 1


def test_rounded_hit_rejects_a_non_number(restore):
    assert stats._rounded_hit("plain text", "1.2.3") is False


def test_literal_source_walks_each_record(tmp_path, restore):
    data = tmp_path / "data"
    data.mkdir()
    root = tmp_path / "root"
    (root / "bbdir").mkdir(parents=True)
    (root / "scripts").mkdir()
    (root / "scripts" / "h2h.py").write_text("seen 8.5\n")
    stats.DATA = data
    stats.ROOT = root
    neigh = data / "neigh-head2head.json"
    neigh.write_text('{"value": 4.5}\n')
    found = stats.literal_source("4.5")
    assert found == ("paper/tex/figures/data/neigh-head2head.json", 1, stats.COMMITS["neigh"])

    log = root / "bbdir" / "campaign-arm.log"
    log.write_text("skip 1\nkept 6.25\n")
    data.joinpath("neigh-head2head.json").unlink()
    found = stats.literal_source("6.25")
    assert found == ("bbdir/campaign-arm.log", 2, stats.COMMITS["measure"])

    log.write_text("nothing here\n")
    found = stats.literal_source("8.5")
    assert found[0] == "scripts/h2h.py"
    assert found[2] == stats.COMMITS["measure"]

    (root / "scripts" / "h2h.py").write_text("no number\n")
    found = stats.literal_source("8.5")
    assert found == ("paper/tex/figures/stats.py", 1, stats.COMMITS["measure"])


def test_emit_numbers_without_the_optional_sidecars(tmp_path, restore):
    root = tmp_path / "root"
    tex = root / "paper" / "tex"
    (tex / "figures" / "data").mkdir(parents=True)
    (tex / "rg_main.tex").write_text("a \\dsv{plain} and \\num{9}\n")
    (tex / "rg_supplement.tex").write_text("b \\qty{1.25}{ms}\n")
    (tex / "figures" / "data" / "tab-scale.tex").write_text("\\num{3}\n")
    data = tmp_path / "data"
    data.mkdir()
    (data / "ring-head2head.json").write_text('{"n": 9, "other": 1.25}\n')
    stats.ROOT = root
    stats.DATA = data
    stats.NUMBER_LOG.clear()
    stats.NUMBER_LOG.append(("plain", "paper/tex/figures/stats.py", 4, "9"))
    dest = tmp_path / "NUMBERS.tsv"
    count = stats.emit_numbers(dest)
    text = dest.read_text()
    assert text.startswith("number\ttex_file\ttex_line\tlog_file\tlog_line\tcode_commit\ton_main\n")
    assert "\t9\t" in text or text.splitlines()[1].startswith("9\t")
    assert count == len(text.splitlines()) - 1
    assert "ring-head2head.json" in text


def test_versions_reject_a_short_ladder(tmp_path, restore):
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="ladder records differ"):
        stats.versions()


def test_versions_reject_the_wrong_call_count(tmp_path, restore):
    stats.DATA = tmp_path
    stats.DRAWS = 4
    for n in stats.SET["ladder"]:
        calls = stats.SET["v1"].get(n, (0, int(stats.SET["v@repeats"])))[1]
        _calls(tmp_path / ("dseams2-ic-%d.txt" % n), int(stats.SET["v@repeats"]), n)
        _err(tmp_path / ("dseams2-ic-%d.err" % n), 2 * n)
        if n in stats.SET["v1"]:
            _calls(tmp_path / ("dseams1-ic-%d.txt" % n), 1, n)
            _err(tmp_path / ("dseams1-ic-%d.err" % n), 2 * n)
    with pytest.raises(SystemExit, match="calls"):
        stats.versions()


def test_versions_reject_the_wrong_six_ring_count(tmp_path, restore):
    stats.DATA = tmp_path
    stats.DRAWS = 4
    _version_tree(tmp_path, bad_six=("two", 216))
    with pytest.raises(SystemExit, match="six-rings on 216"):
        stats.versions()


def test_versions_reject_a_short_mw_sample(tmp_path, restore):
    stats.DATA = tmp_path
    stats.DRAWS = 4
    _version_tree(tmp_path, mw_rows=2)
    with pytest.raises(SystemExit, match="mW crystal"):
        stats.versions()


def test_versions_reject_the_wrong_mw_census(tmp_path, restore):
    stats.DATA = tmp_path
    stats.DRAWS = 4
    _version_tree(tmp_path, mw_atoms=8, mw_six=1)
    with pytest.raises(SystemExit, match="mW molecules"):
        stats.versions()


def test_churn_rejects_a_different_protocol(tmp_path, restore):
    text = (DATA_SRC / "churn-sweep.txt").read_text()
    (tmp_path / "churn-sweep.txt").write_text(text.replace("frames 147->148", "frames 1->2"))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="sweep record differs"):
        stats.churn()


def test_churn_rejects_a_call_that_misses_the_row(tmp_path, restore):
    lines = (DATA_SRC / "churn-sweep.txt").read_text().splitlines()
    for index, line in enumerate(lines):
        if line.startswith("calls "):
            head = line.split()[:3]
            lines[index] = " ".join(head + ["1"])
            break
    (tmp_path / "churn-sweep.txt").write_text("\n".join(lines) + "\n")
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="churn calls disagree"):
        stats.churn()


def _h2h_copy(tmp_path):
    data = tmp_path / "data"
    shutil.copytree(DATA_SRC / "h2h", data / "h2h")
    stats.DATA = data
    stats.H2H = data / "h2h"
    return data


def test_h2h_covers_the_gaps_in_the_records(tmp_path, restore):
    data = _h2h_copy(tmp_path)
    path = data / "h2h" / "comparisons.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    fields = list(rows[0])
    kept = []
    for row in rows:
        if row["system"] == "ic" and row["scale"] == "1e3" and row["other"] == "genice-rings":
            continue
        if row["other"].startswith("matscipy-rings-"):
            row["equivalent_workload"] = "False"
        kept.append(row)
    extra = dict(kept[0])
    extra["other"] = "extra-code"
    extra["equivalent_workload"] = "True"
    extra["other_over_dseams"] = "4"
    kept.append(extra)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)
    timings = data / "h2h" / "timings.tsv"
    timing_rows = list(csv.DictReader(timings.open(), delimiter="\t"))
    timing_fields = list(timing_rows[0])
    timing_rows = [
        row for row in timing_rows
        if not (row["system"] == "liquid" and row["scale"] == "1e6")
    ]
    with timings.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, timing_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(timing_rows)
    stats.macros.clear()
    stats.h2h()
    text = _macro_text()
    assert "dsv@hl@cycless@ih\\endcsname" in text
    assert "dsv@hl@cycless\\endcsname" not in text
    assert "dsv@h@range@matscipy@min\\endcsname" not in text
    ladder = (data / "tab-ladder.tex").read_text()
    assert "884736" not in ladder
    assert "Ice~Ic" in ladder


def test_worked_ignores_a_line_that_is_not_a_count(tmp_path, restore):
    data = tmp_path / "data"
    data.mkdir()
    shutil.copy(DATA_SRC / "worked_example.tsv", data / "worked_example.tsv")
    ideal = (DATA_SRC / "ideal-cages.txt").read_text() + "\nnote extra\n"
    (data / "ideal-cages.txt").write_text(ideal)
    update = (DATA_SRC / "worked_update.tsv").read_text() + "note\tkept\n"
    (data / "worked_update.tsv").write_text(update)
    stats.DATA = data
    stats.macros.clear()
    stats.worked()
    text = _macro_text()
    assert "dsv@w@ice-Ic@sites\\endcsname{\\num{1000}}" in text
    assert "dsv@u@sites\\endcsname{\\num{72}}" in text
    assert "dsv@u@note\\endcsname" not in text


def test_nucleation_keeps_a_comment_out_of_the_frame_count(tmp_path, restore):
    data = tmp_path / "data"
    data.mkdir()
    for name in (
        "nucleation-cutoff-carry.txt",
        "nucleation-seeded-carry.txt",
        "nucleation-host.txt",
        "nucleation-timings.txt",
    ):
        shutil.copy(DATA_SRC / name, data / name)
    carry = data / "nucleation-cutoff-carry.txt"
    carry.write_text("# comment\n" + carry.read_text())
    stats.DATA = data
    stats.DRAWS = 3
    stats.macros.clear()
    stats.nucleation()
    assert "dsv@n@frames\\endcsname{148}" in _macro_text()


def test_nucleation_rejects_a_short_trajectory(tmp_path, restore):
    data = tmp_path / "data"
    data.mkdir()
    (data / "nucleation-cutoff-carry.txt").write_text("nop 4 graph cutoff hexagonal 1 cubic 1 both 0 water 2\n")
    (data / "nucleation-seeded-carry.txt").write_text("nop 4 graph seeded hexagonal 1 cubic 1 both 0 water 2\n")
    stats.DATA = data
    with pytest.raises(SystemExit, match="nucleation runs read 1 frames"):
        stats.nucleation()


def test_nucleation_rejects_a_short_timing_sample(tmp_path, restore):
    data = tmp_path / "data"
    data.mkdir()
    for name in ("nucleation-cutoff-carry.txt", "nucleation-seeded-carry.txt", "nucleation-host.txt"):
        shutil.copy(DATA_SRC / name, data / name)
    (data / "nucleation-timings.txt").write_text("cutoff carry 1 1.5\n")
    stats.DATA = data
    with pytest.raises(SystemExit, match="timed runs"):
        stats.nucleation()


def test_ovito_rejects_a_short_sample(tmp_path, restore):
    data = tmp_path / "data"
    h2h = data / "h2h"
    h2h.mkdir(parents=True)
    text = (DATA_SRC / "h2h" / "ovito-rings.tsv").read_text().splitlines()
    head, row = text[0], text[1].split("\t")
    row[-1] = "1.0,2.0"
    (h2h / "ovito-rings.tsv").write_text(head + "\n" + "\t".join(row) + "\n")
    stats.DATA = data
    stats.H2H = h2h
    with pytest.raises(SystemExit, match="timed calls"):
        stats.ovito_rings()


def _fake_rings():
    rings = {}
    for system in stats.H2H_SYSTEMS:
        for scale in ("1e3", "1e4", "1e5", "1e6"):
            rings[(system, scale, "dseams-rings")] = [8, 8, 8, 8, 8, 8]
            rings[(system, scale, "ovito-rings")] = [24, 24, 24, 24, 24, 24]
    return rings


def test_h2h_table_prints_a_ratio_and_a_gap(tmp_path, restore):
    data = _h2h_copy(tmp_path)
    stats.DRAWS = 4
    comps = data / "h2h" / "comparisons.tsv"
    rows = list(csv.DictReader(comps.open(), delimiter="\t"))
    fields = list(rows[0])
    rows = [
        row for row in rows
        if not (row["system"] == "liquid" and row["scale"] == "1e5" and row["other"] == "ovito-chill")
    ]
    added = dict(rows[0])
    added.update({
        "system": "ic",
        "scale": "1e5",
        "other": "genice-rings",
        "other_over_dseams": "12.5",
        "ratio_ci95_low": "12",
        "ratio_ci95_high": "13",
    })
    rows.append(added)
    with comps.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    errors = data / "h2h" / "errors.jsonl"
    kept = []
    for line in errors.read_text().splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("system") == "liquid" and row.get("scale") == "1e5" and row.get("method") == "ovito-chill":
            continue
        kept.append(line)
    errors.write_text("\n".join(kept) + "\n")
    timings = data / "h2h" / "timings.tsv"
    timing_rows = list(csv.DictReader(timings.open(), delimiter="\t"))
    timing_fields = list(timing_rows[0])
    timing_rows = [
        row for row in timing_rows
        if not (row["system"] == "liquid" and row["scale"] == "1e5" and row["method"] == "dseams-chill")
    ]
    with timings.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, timing_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(timing_rows)
    stats.macros.clear()
    stats.h2h_table(_fake_rings())
    table = (data / "tab-h2h.tex").read_text()
    assert "not rerun" in table
    assert "12.5" in table
    assert "rings & d-SEAMS" in table
    assert "\\num{8.00}" in table


def test_ranges_reject_an_unexpected_matched_set(tmp_path, restore):
    data = _h2h_copy(tmp_path)
    path = data / "h2h" / "comparisons.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    for row in rows:
        if row["other"] == "ovito-chill":
            row["equivalent_workload"] = "False"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(SystemExit, match="matched codes differ"):
        stats.ranges(_fake_rings())


def test_ranges_leave_a_missing_code_as_not_rerun(tmp_path, restore):
    data = _h2h_copy(tmp_path)
    path = data / "h2h" / "comparisons.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    rows = [row for row in rows if row["other"] != "ovito-ptm"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    stats.macros.clear()
    stats.ranges(_fake_rings())
    table = (data / "tab-ranges.tex").read_text()
    assert "not rerun" in table
    assert "dsv@r@matched\\endcsname" in _macro_text()


def test_ranges_reject_a_matched_code_that_is_not_slower(tmp_path, restore):
    data = _h2h_copy(tmp_path)
    path = data / "h2h" / "comparisons.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    rows[0]["other_over_dseams"] = "0.5"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(SystemExit, match="not slower"):
        stats.ranges(_fake_rings())


def test_thresholds_reject_a_missing_default(tmp_path, restore):
    source = tmp_path / "src"
    source.mkdir()
    (source / "cage_affiliation.cpp").write_text("no default\n")
    stats.SOURCE = source
    stats.THRESHOLD_COPIES = ()
    with pytest.raises(SystemExit, match="0 matches"):
        stats.thresholds()


def test_thresholds_reject_two_defaults(tmp_path, restore):
    source = tmp_path / "src"
    source.mkdir()
    (source / "cage_affiliation.cpp").write_text(
        "double batchFraction = 0.25;\ndouble batchFraction = 0.25;\n"
    )
    stats.SOURCE = source
    with pytest.raises(SystemExit, match="2 matches"):
        stats.thresholds()


def test_a_script_must_copy_the_src_default(tmp_path, restore):
    (tmp_path / "copy.cpp").write_text("double batchFraction = 0.25;\n")
    stats.ROOT = tmp_path
    stats.THRESHOLD_COPIES = (("cagerows", "copy.cpp", r"double batchFraction = ([0-9.]+);"),)
    stats.macros.clear()
    stats.thresholds()
    assert "dsv@thr@cagerows\\endcsname{0.25}" in _macro_text()
    (tmp_path / "copy.cpp").write_text("double batchFraction = 0.5;\n")
    with pytest.raises(SystemExit, match="another value"):
        stats.thresholds()


def test_read_setting_requires_one_shared_value(tmp_path, restore):
    (tmp_path / "a.txt").write_text("note\n")
    stats.ROOT = tmp_path
    with pytest.raises(SystemExit, match="no match"):
        stats.read_setting("cutoff", (("a.txt", r"^CUT = ([0-9.]+)$"),))
    (tmp_path / "a.txt").write_text("CUT = 3.5\n")
    (tmp_path / "b.txt").write_text("CUT = 4\n")
    with pytest.raises(SystemExit, match="disagree"):
        stats.read_setting("cutoff", (
            ("a.txt", r"^CUT = ([0-9.]+)$"),
            ("b.txt", r"^CUT = ([0-9.]+)$"),
        ))
    assert stats.read_setting("cutoff", (("a.txt", r"^CUT = ([0-9.]+)$"),)) == decimal.Decimal("3.5")
    assert stats.source_line("a.txt", r"^ABSENT") == 1
    assert stats.source_line("a.txt", r"^CUT = ") == 1


def test_settings_reject_a_ring_column_that_moved(restore):
    stats.MIN_RING = 8
    with pytest.raises(SystemExit, match="ring columns"):
        stats.settings()


def test_settings_reject_a_different_prism_list(restore):
    stats.PRISM_SIZES = (4, 5)
    with pytest.raises(SystemExit, match="prism columns"):
        stats.settings()


def _settings_on(tmp_path, edit):
    root = tmp_path / "root"
    _plant(root)
    edit(root)
    stats.ROOT = root
    stats.macros.clear()
    stats.settings()


def test_settings_reject_a_campaign_that_dropped_the_ring_bound(tmp_path, restore):
    def edit(root):
        path = root / "scripts" / "h2h.py"
        path.write_text(path.read_text().replace("ringNetwork(idx, 6)", "ringNetwork(idx, 4)"))

    with pytest.raises(SystemExit, match="no longer counts the rings"):
        _settings_on(tmp_path, edit)


def test_settings_reject_a_broken_ladder(tmp_path, restore):
    def edit(root):
        path = root / "paper" / "tex" / "repro" / "versions.sh"
        path.write_text(path.read_text().replace('CELLS="3 4 5 6 7 8"', 'CELLS="3 5 8"'))

    with pytest.raises(SystemExit, match="consecutive ladder cells"):
        _settings_on(tmp_path, edit)


def test_settings_reject_uneven_dseams1_calls(tmp_path, restore):
    def edit(root):
        path = root / "paper" / "tex" / "repro" / "versions.sh"
        path.write_text(path.read_text().replace("1728:1:6", "1728:2:9"))

    with pytest.raises(SystemExit, match="different calls"):
        _settings_on(tmp_path, edit)


def test_settings_reject_a_split_ovito_plan(tmp_path, restore):
    def edit(root):
        path = root / "paper" / "tex" / "figures" / "stats.py"
        path.write_text(path.read_text().replace(
            'PLAN = {"1e3": (5, 30), "1e4": (5, 30), "1e5": (5, 30), "1e6": (1, 10)}',
            'PLAN = {"1e3": (5, 30), "1e4": (1, 10), "1e6": (1, 10)}',
        ))

    with pytest.raises(SystemExit, match="varies its calls"):
        _settings_on(tmp_path, edit)


def test_hardware_rejects_a_different_processor(tmp_path, restore):
    data = tmp_path / "data"
    h2h = data / "h2h"
    h2h.mkdir(parents=True)
    shutil.copy(DATA_SRC / "host.txt", data / "host.txt")
    churn = (DATA_SRC / "churn-host.txt").read_text().replace("Intel Xeon Processor", "Other CPU")
    (data / "churn-host.txt").write_text(churn)
    shutil.copy(DATA_SRC / "h2h" / "environment.json", h2h / "environment.json")
    stats.DATA = data
    stats.H2H = h2h
    with pytest.raises(SystemExit, match="machine records differ"):
        stats.hardware()


def test_marks_reject_a_spread_outside_the_mark(tmp_path, restore):
    data = tmp_path / "data"
    data.mkdir()
    (data / "plot-versions-x.dat").write_text(
        "sites median lo hi q1 q3\n10 1 0.1 1 1 1\n"
    )
    stats.DATA = data
    with pytest.raises(SystemExit, match="leaves the marks"):
        stats.marks()


def test_a_figure_without_records_is_rejected(restore, monkeypatch):
    real = Path.read_text

    def read_text(self, *args, **kwargs):
        if self.name == "fig-versions.tex":
            return "no records\n"
        return real(self, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", read_text)
    with pytest.raises(SystemExit, match="plots no records"):
        stats.figure_calls()


def _load_json(name, dest):
    record = json.loads((DATA_SRC / name).read_text())
    dest.write_text(json.dumps(record))
    return record


def test_panels_are_skipped_when_the_record_is_absent(tmp_path, restore):
    stats.DATA = tmp_path
    assert stats.neigh_panel() is None
    assert stats.ring_panel() is None


def test_neigh_panel_rejects_a_second_thread(tmp_path, restore):
    record = json.loads((DATA_SRC / "neigh-head2head.json").read_text())
    record["thread_policy"]["ovito_threads"] = 4
    (tmp_path / "neigh-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="ovito_threads is 4"):
        stats.neigh_panel()


def test_neigh_panel_rejects_extra_dseams_threads(tmp_path, restore):
    record = json.loads((DATA_SRC / "neigh-head2head.json").read_text())
    record["thread_policy"]["dseams_threads_reported"] = [1, 2]
    (tmp_path / "neigh-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="not one d-SEAMS thread"):
        stats.neigh_panel()


def test_neigh_panel_rejects_the_wrong_repeat_count(tmp_path, restore):
    record = json.loads((DATA_SRC / "neigh-head2head.json").read_text())
    record["protocol"]["repeats"] = 11
    (tmp_path / "neigh-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="not 12 timed calls"):
        stats.neigh_panel()


def test_neigh_panel_skips_a_stage_that_is_not_a_table(tmp_path, restore):
    record = json.loads((DATA_SRC / "neigh-head2head.json").read_text())
    record["datasets"]["orth-4k"]["knn_stage"] = None
    (tmp_path / "neigh-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    stats.macros.clear()
    assert stats.neigh_panel() is True
    assert "dsv@neigh@skew@count\\endcsname" in _macro_text()
    assert (tmp_path / "tab-neigh-perf.tex").is_file()


def test_neigh_panel_rejects_a_short_sample(tmp_path, restore):
    record = json.loads((DATA_SRC / "neigh-head2head.json").read_text())
    stage = record["datasets"]["orth-4k"]["cutoff_stage"]
    tool = next(iter(stage))
    stage[tool]["samples_ms"] = [1.0]
    (tmp_path / "neigh-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="does not have 12 samples"):
        stats.neigh_panel()


def test_neigh_panel_rejects_disagreeing_skew_counts(tmp_path, restore):
    record = json.loads((DATA_SRC / "neigh-head2head.json").read_text())
    record["datasets"]["tri-skew-512"]["knn_stage"]["dseams_omp1"]["count"] = 1
    (tmp_path / "neigh-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="skew neighbour counts differ"):
        stats.neigh_panel()


def test_ring_panel_rejects_a_second_thread(tmp_path, restore):
    record = json.loads((DATA_SRC / "ring-head2head.json").read_text())
    record["thread_policy"]["ovito_threads"] = 2
    (tmp_path / "ring-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="not one thread"):
        stats.ring_panel()


def test_ring_panel_rejects_the_wrong_protocol(tmp_path, restore):
    record = json.loads((DATA_SRC / "ring-head2head.json").read_text())
    record["protocol"]["warmup_runs"] = 1
    (tmp_path / "ring-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="not 12 timed calls"):
        stats.ring_panel()


def test_ring_panel_rejects_disagreeing_counts(tmp_path, restore):
    record = json.loads((DATA_SRC / "ring-head2head.json").read_text())
    record["datasets"]["tri-136"]["ring_stage_ms"]["ovito"]["count"] = 1
    (tmp_path / "ring-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="ring counts differ"):
        stats.ring_panel()


def test_ring_panel_rejects_a_short_sample(tmp_path, restore):
    record = json.loads((DATA_SRC / "ring-head2head.json").read_text())
    record["datasets"]["tri-136"]["ring_stage_ms"]["ovito"]["samples_ms"] = [1.0]
    (tmp_path / "ring-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    with pytest.raises(SystemExit, match="does not have 12 samples"):
        stats.ring_panel()


def test_ring_panel_omits_setup_and_bonds_when_they_were_not_timed(tmp_path, restore):
    record = json.loads((DATA_SRC / "ring-head2head.json").read_text())
    record["datasets"]["tri-136"].pop("setup_ms")
    record["datasets"]["orth-4k"].pop("bonds_and_rings_ms")
    (tmp_path / "ring-head2head.json").write_text(json.dumps(record))
    stats.DATA = tmp_path
    quoted = stats.ring_panel()
    table = (tmp_path / "tab-ring-perf.tex").read_text()
    assert table.count("setup, cold") == 2
    assert table.count("bonds and rings") == 2
    assert 136 in quoted
    assert (tmp_path / "plot-serial-rings.dat").read_text().splitlines()[0] == "sites dseams ovito"


def _prepare_main(tmp_path):
    data = tmp_path / "data"
    shutil.copytree(DATA_SRC, data)
    stats.DATA = data
    stats.H2H = data / "h2h"
    stats.DRAWS = 5
    stats.macros.clear()
    stats.NUMBER_LOG.clear()
    stats.plot_calls.clear()
    return data


def test_a_macro_defined_twice_stops_the_driver(tmp_path, restore):
    _prepare_main(tmp_path)
    real = stats.define

    def wrapped(name, body):
        real(name, body)
        if name == "set@cutoff" and wrapped.once:
            wrapped.once = False
            real(name, "0")

    wrapped.once = True
    stats.define = wrapped
    with pytest.raises(SystemExit, match="defined twice"):
        stats.main()


def test_the_campaigns_must_count_the_same_matched_pairs(tmp_path, restore):
    data = _prepare_main(tmp_path)
    path = data / "h2h" / "comparisons.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    flipped = False
    for row in rows:
        if row["equivalent_workload"] == "True" and row["other"] == "ovito-chill":
            row["equivalent_workload"] = "False"
            flipped = True
            break
    assert flipped
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(SystemExit, match="different matched pairs"):
        stats.main()


def test_main_writes_the_published_numbers(tmp_path, restore, capsys):
    data = tmp_path / "data"
    shutil.copytree(DATA_SRC, data)
    argv = sys.argv[:]
    sys.argv = ["stats.py", str(data)]
    try:
        runpy.run_path(str(REPO / "paper" / "tex" / "figures" / "stats.py"), run_name="__main__")
    finally:
        sys.argv = argv
    out = capsys.readouterr().out
    assert "4897 macros" in out
    assert "numbers 1092" in out
    assert "x86_64 | 8 | Intel Xeon Processor" in out
    macros = (data / "macros.tex").read_text()
    assert "dsv@set@q6\\endcsname{0.558}" in macros
    assert "dsv@set@repeats\\endcsname{6}" in macros
    assert "dsv@v@216@ring@one\\endcsname{4.61}" in macros
    assert "dsv@hl@matched\\endcsname{79}" in macros
    versions = (data / "tab-versions.tex").read_text()
    assert "\\num{2744}" in versions
    assert "not run" in versions
    numbers = REPO / "bbdir" / "NUMBERS.tsv"
    assert numbers.read_text().startswith("number\t")
    assert sum(1 for _ in numbers.open()) == 1093
