# Copyright (c) d-SEAMS developers
"""Cover the acceptance helpers and the error branches around them."""

import os
import shutil
import stat
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_build = os.environ.get("BUILD", "")
if _build:
    os.environ.setdefault("SEAMS_LIB", _build + "/seams/src/libyodaLib.so")
    os.environ.setdefault("SEAMS_CLI", _build + "/seams/seams")
    os.environ.setdefault("SEAMS_LAMMPS_PLUGIN", _build + "/plugin/libdseams_lammps.so")
    os.environ.setdefault("DSEAMS_MD_PREFIX", _build + "/md-prefix")
os.environ.setdefault("DSEAMS_LIVE_WORK", "var/run/dseams-live-cover")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import engine_e2e  # noqa: E402
import live  # noqa: E402


def expect_error(fn):
    try:
        fn()
    except (RuntimeError, SystemExit):
        return True
    return False


def check(cond):
    if not cond:
        raise SystemExit("check failed")


def write_dump(path, step=1):
    path.write_text(
        "ITEM: TIMESTEP\n%d\nITEM: NUMBER OF ATOMS\n1\n"
        "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
        "ITEM: ATOMS id type x y z\n1 1 0 0 0\n" % step
    )


def real_lines(dump):
    rings = [ln for ln in live.cli_text("rings", dump).splitlines() if ln.startswith("rings ")][0]
    chill = [ln for ln in live.cli_text("chill-plus", dump).splitlines() if ln.startswith("nop ")][0]
    cage = [
        ln
        for ln in live.cli_text("cages", dump, ["--graph", "cutoff"]).splitlines()
        if ln.startswith("nop ")
    ][0]
    return rings, chill, cage


def main():
    work = Path(os.environ["DSEAMS_LIVE_WORK"])
    work.mkdir(parents=True, exist_ok=True)
    live.WORK = work
    try:
        check(False)
    except SystemExit:
        pass
    check(True)
    quiet_call = expect_error(lambda: None)
    check(quiet_call is False)
    check(expect_error(lambda: live.run(["/bin/false"], work)))
    folder = work / "unit"
    folder.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(folder)
    folder.mkdir()
    dump = folder / "one.lammpstrj"
    write_dump(dump)
    record = folder / "live.txt"

    quiet = live.run(["/bin/false"], folder, check=False)
    check(quiet.returncode != 0)

    info = live.host_block()
    check(bool(info["head"]) and bool(info["cpu"]))
    digest = live.sha256(dump)
    check(len(digest) == 64)

    blocks = live.census_blocks("note\nstep 1\nrings 3:0\nchill cubic:1\ncage hc:0\n")
    check(blocks[0]["step"] == 1 and blocks[0]["chill"]["cubic"] == 1)

    live.lammps_input(folder / "in-dump", 1, dump, record, "plugin.so", write_dump=True)
    live.lammps_input(folder / "in-nodump", 1, dump, record, "plugin.so", write_dump=False)
    live.lammps_plain(folder / "plain-nodump", 1, every=0)
    live.lammps_plain(folder / "plain-dump", 1, every=1, native=dump)

    short = folder / "short.dcd"
    short.write_bytes(b"1234")
    check(expect_error(lambda: live.native_frame_count(short)))
    check(expect_error(lambda: live.dcd_schedule(short)))
    import struct
    shifted = folder / "shifted.dcd"
    shifted.write_bytes(struct.pack("<i4siiii", 84, b"CORD", 2, 10, 5, 15) + b"\0" * 64)
    record.write_text("step 5\nrings 3:0\nstep 10\nrings 3:0\n")
    shifted_error = False
    try:
        live.match_native_cli(record, shifted, type_filter=8)
    except RuntimeError as exc:
        shifted_error = "DCD step list differs" in str(exc)
    check(shifted_error)
    matched = folder / "matched.dcd"
    matched.write_bytes(struct.pack("<i4siiii", 84, b"CORD", 2, 5, 5, 10) + b"\0" * 64)
    check(live.dcd_schedule(matched) == [5, 10])
    one = folder / "one.dcd"
    one.write_bytes(struct.pack("<i4siiii", 84, b"CORD", 1, 5, 5, 5) + b"\0" * 64)
    record.write_text("step 1\nrings 3:0\n")
    one_error = False
    try:
        live.match_native_cli(record, one, type_filter=8)
    except RuntimeError as exc:
        one_error = "DCD step list differs" in str(exc)
    check(one_error)
    bad_stride = folder / "bad-stride.dcd"
    bad_stride.write_bytes(struct.pack("<i4siiii", 84, b"CORD", 2, 5, 0, 10) + b"\0" * 64)
    check(expect_error(lambda: live.dcd_schedule(bad_stride)))
    check(live.native_frame_count(dump) == 1)

    bad_trr = folder / "bad.trr"
    bad_trr.write_bytes(b"\x00" * 84)
    check(expect_error(lambda: live.trr_frames(bad_trr)))
    import struct
    header = bytearray(b"\x00" * 84)
    struct.pack_into(">i", header, 0, 1993)
    struct.pack_into(">i", header, 32, 36)
    struct.pack_into(">i", header, 52, 0)
    struct.pack_into(">ii", header, 64, 0, 1)
    (folder / "trunc.trr").write_bytes(bytes(header) + b"\x00" * 40)
    check(expect_error(lambda: live.trr_frames(folder / "trunc.trr")))

    check(live.gmx_wall("Time: 1.5 2.5\n") == 2.5)
    check(expect_error(lambda: live.gmx_wall("no time here\n")))

    check(expect_error(lambda: live.ladder_values("1\n2\n", "")))
    samples = live.ladder_values("\n".join(["0.1"] * 7) + "\n", "")
    check(len(samples) == 7)

    check(expect_error(lambda: live.require_equal(1, 2, "mismatch")))
    live.require_equal(3, 3, "same")
    check(expect_error(lambda: live.require_native_frames(1, 2, 2)))
    check(expect_error(lambda: live.require_native_frames(1, 1, 2)))
    live.require_native_frames(4, 4, 4)

    rings, chill, cage = real_lines(dump)
    chill_bits = []
    chill_words = chill.split()
    for index in range(2, len(chill_words), 2):
        chill_bits.append("%s:%s" % (chill_words[index], chill_words[index + 1]))
    cage_words = cage.split()
    cage_map = {cage_words[index]: cage_words[index + 1] for index in range(4, len(cage_words), 2)}
    cage_bits = "hc:%s ddc:%s both:%s water:%s" % (
        cage_map.get("hexagonal", "0"),
        cage_map.get("cubic", "0"),
        cage_map.get("both", "0"),
        cage_map.get("water", "0"),
    )
    good = "step 1\n%s\nchill %s\ncage %s\n" % (rings, " ".join(chill_bits), cage_bits)
    record.write_text("step 1\nstep 1\n" + rings + "\n")
    check(expect_error(lambda: live.match_lammps_cli(record, dump)))
    write_dump(folder / "dup.lammpstrj", 1)
    (folder / "dup.lammpstrj").write_text(dump.read_text() + dump.read_text())
    record.write_text("step 1\n" + rings + "\nchill cubic:0\ncage hc:0 ddc:0 both:0 water:0\n")
    check(expect_error(lambda: live.match_lammps_cli(record, folder / "dup.lammpstrj")))
    record.write_text("step 9\n" + rings + "\n")
    check(expect_error(lambda: live.match_lammps_cli(record, dump)))

    saved = live.SEAMS
    stub = folder / "stub-seams"
    stub.write_text("#!/bin/sh\nexit 0\n")
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    live.SEAMS = str(stub.resolve())
    record.write_text("step 1\n" + rings + "\n")
    check(expect_error(lambda: live.match_lammps_cli(record, dump)))
    check(expect_error(lambda: live.match_native_cli(record, dump)))
    live.SEAMS = saved

    record.write_text("step 1\nrings not-the-census\nchill cubic:0\ncage hc:0 ddc:0 both:0 water:0\n")
    check(expect_error(lambda: live.match_lammps_cli(record, dump)))
    record.write_text("step 1\n" + rings + " extra\nchill cubic:9\ncage hc:0 ddc:0 both:0 water:0\n")
    check(expect_error(lambda: live.match_lammps_cli(record, dump)))
    record.write_text("step 1\n" + rings + "\nchill " + " ".join(chill_bits) + "\ncage hc:9 ddc:0 both:0 water:0\n")
    check(expect_error(lambda: live.match_lammps_cli(record, dump)))
    record.write_text("step 1\n" + rings + "\nchill " + " ".join(chill_bits) + "\ncage hc:0 ddc:0 both:0 water:9\n")
    check(expect_error(lambda: live.match_lammps_cli(record, dump)))
    record.write_text("step 1\n" + rings + "\nchill " + " ".join(chill_bits) + "\ncage " + cage_bits + "\n")
    record.write_text("step 1\n" + rings + "\nchill " + " ".join(chill_bits) + "\ncage hc:9 ddc:0 both:0 water:0\n")
    check(expect_error(lambda: live.match_native_cli(record, dump)))
    record.write_text("step 1\nrings 3:9 4:0 5:0 6:0 7:0 8:0\nchill " + " ".join(chill_bits) + "\ncage " + cage_bits + "\n")
    check(expect_error(lambda: live.match_native_cli(record, dump)))

    parsed = live.parse_dump(dump)
    check(len(parsed) == 1 and parsed[0][0] == 1)
    (folder / "empty.txt").write_text("")
    check(live.parse_dump(folder / "empty.txt") == [])

    check(expect_error(lambda: live.module_entry(["live.py"])))
    live.module_entry([sys.executable, "--ladder-worker", str(dump)])

    empty = folder / "empty.lammpstrj"
    empty.write_text("")
    check(expect_error(lambda: live.thread_ladder(empty, 1)))
    live.thread_ladder(dump.resolve(), 1)

    made = live.diamond_positions(1)
    lattice = work / "lattice" / "1" / "one.lammpstrj"
    lattice.unlink()
    made = live.diamond_positions(1)
    check(len(made) >= 1)

    check(expect_error(lambda: engine_e2e.need("DSEAMS_NOT_SET")))
    saved_mpi = os.environ.pop("DSEAMS_MPIRUN", None)
    check(expect_error(lambda: engine_e2e.lammps_case("unused", "plugin.so", work / "rank-gap", 2)))
    if saved_mpi is not None:
        os.environ["DSEAMS_MPIRUN"] = saved_mpi

    stub_lmp = folder / "stub-lmp"
    stub_lmp.write_text(
        "#!/bin/sh\n"
        "cat > native.lammpstrj << 'END'\n"
        + dump.read_text()
        + "END\n"
        "cat > live.txt << 'END'\n"
        + good
        + "END\n"
    )
    stub_lmp.chmod(stub_lmp.stat().st_mode | stat.S_IEXEC)
    check(expect_error(lambda: engine_e2e.lammps_case(str(stub_lmp.resolve()), "plugin.so", work / "short-run", 1)))

    for name in ("lammps-plugin-1", "lammps-plugin-2", "lammps-package", "e2e-gromacs"):
        (work / name).mkdir(parents=True, exist_ok=True)
    check(live.SEAMS == os.environ["SEAMS_CLI"])
    from dseams_insitu.engine import Engine
    with Engine(os.environ["SEAMS_LIB"], str(folder / "path-record.txt")) as eng:
        eng.observe(0, [[0.0, 0.0, 0.0], [2.8, 0.0, 0.0], [0.0, 2.8, 0.0], [1.4, 1.4, 2.0]],
                    [0.0, 0.0, 0.0, 20.0, 20.0, 20.0, 0.0, 0.0, 0.0], None, [1, 1, 1, 1])
    os.environ["DSEAMS_E2E_WORK"] = str(work)
    engine_e2e.main()
    check(live.type7_median([1.0, 3.0]) == 2.0)
    gro = folder / "box.gro"
    gro.write_text("title\n1\n    1SOL     OW    1   0.1   0.1   0.1\n   2.0   2.0   2.0\n")
    check(live.gromacs_box(gro)[3] == 20.0)
    check(live.free_port() > 0)
    check(len(live.ice_positions()) == 8)
    saved_exe = os.environ.get("SEAMS_INSITU_EXE")
    saved_imd = os.environ.get("SEAMS_GRO_IMD_EXE")
    os.environ["SEAMS_INSITU_EXE"] = ""
    os.environ["SEAMS_GRO_IMD_EXE"] = ""
    check(live.acceptance() is None)
    if saved_exe:
        os.environ["SEAMS_INSITU_EXE"] = saved_exe
    else:
        os.environ.pop("SEAMS_INSITU_EXE", None)
    if saved_imd:
        os.environ["SEAMS_GRO_IMD_EXE"] = saved_imd
    else:
        os.environ.pop("SEAMS_GRO_IMD_EXE", None)
    for name in (
        "run_lammps",
        "ase_identity",
        "run_ipi",
        "gromacs_trr_identity",
        "lammps_mpi_identity",
        "gromacs_mpi_identity",
        "openmm_identity",
    ):
        result = getattr(live, name)()
        if name == "run_lammps":
            frames = result.get("frames") if isinstance(result, dict) else None
            if isinstance(frames, bool) or not isinstance(frames, int) or frames < 1:
                raise SystemExit(name + " returned %r" % (result,))
        elif name in ("lammps_mpi_identity", "gromacs_mpi_identity"):
            if not (isinstance(result, tuple) and len(result) == 2 and result[0] is True and isinstance(result[1], str)):
                raise SystemExit(name + " returned %r" % (result,))
        elif isinstance(result, bool) or not isinstance(result, int) or result < 1:
            raise SystemExit(name + " returned %r" % (result,))
        print(name, "passed", flush=True)
    check(live.ipi_elapsed(False) > 0)
    check(live.ipi_elapsed(True) > 0)
    os.environ.pop("DSEAMS_SCALE_TINY", None)
    check(live.scale_plan()[0][1] == 1000)
    os.environ["DSEAMS_SCALE_TINY"] = "1"
    dup = folder / "dup.txt"
    dup.write_text("step 1\nrings 3:0\nstep 1\nrings 3:0\n")
    check(expect_error(lambda: live.match_lammps_cli(dup, dump)))
    same_a = folder / "same-a.txt"
    same_b = folder / "same-b.txt"
    same_a.write_text("step 1\nrings 3:0\n")
    same_b.write_text(same_a.read_text())
    check(live.same(same_a, same_b) >= 1)
    differ = folder / "differ.txt"
    differ.write_text("step 2\nrings 3:0\n")
    check(expect_error(lambda: live.same(same_a, differ)))
    empty = folder / "empty.txt"
    empty.write_text("no census\n")
    check(expect_error(lambda: live.same(empty, empty)))
    check(expect_error(lambda: live.loop_time("no loop time here\n")))
    blank_xyz = folder / "blank.xyz"
    blank_xyz.write_text("\n1\ncomment\nO 0 0 0\n")
    check(len(live.read_multi_xyz(blank_xyz)) == 1)
    bad_gro = folder / "bad.gro"
    bad_gro.write_text("title\n0\n1.0\n")
    check(expect_error(lambda: live.gromacs_box(bad_gro)))
    check(live.type7_median([4.0]) == 4.0)
    saved_main = live.main
    def boom():
        raise RuntimeError("forced")
    live.main = boom
    os.environ["SEAMS_INSITU_EXE"] = saved_exe or "present"
    os.environ["SEAMS_GRO_IMD_EXE"] = saved_imd or "present"
    check(live.acceptance() is None)
    live.main = saved_main
    check(live.module_entry(["live.py", "--acceptance"]) == 0)
    saved_gmx = os.environ.pop("DSEAMS_GMX", None)
    note = live.gromacs_mpi_identity()
    check(isinstance(note, tuple) and len(note) == 2 and note[0] is False and isinstance(note[1], str) and note[1] != "")
    if saved_gmx is not None:
        os.environ["DSEAMS_GMX"] = saved_gmx
    check(live.loop_time("Loop time of 1.5 on 1 procs for 2 steps\n") == 1.5)
    os.environ["DSEAMS_SCALE_TINY"] = "1"
    live.time_lammps("plain", 2, 1, live.PLUGIN)
    live.time_lammps("live", 2, 1, live.PLUGIN)
    live.main()
    import test_lammps_step0
    test_lammps_step0.run()
    print("cover ok")


if __name__ == "__main__":
    main()
