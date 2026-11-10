# Copyright (c) d-SEAMS developers
"""Compare dumped host neighbour sets with freud and OVITO.

The fix writes the neighbour rows it stored. The compute writes its
per-atom labels. This script reads those text files. It does not
import d-SEAMS. Atomic style has one atom per molecule, so the
per-molecule neighbour set is the per-atom set.
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np

try:
    import freud
except ImportError:
    freud = None

try:
    import ovito
    from ovito.data import CutoffNeighborFinder, DataCollection, Particles, SimulationCell
    from ovito.modifiers import ChillPlusModifier
except ImportError:
    ovito = None
    CutoffNeighborFinder = None
    ChillPlusModifier = None


DSEAMS_NAMES = (
    "cubic",
    "hexagonal",
    "water",
    "interfacial",
    "clathrate",
    "interClathrate",
    "unclassified",
    "reCubic",
    "reHex",
)
OVITO_NAMES = (
    "OTHER",
    "HEXAGONAL",
    "CUBIC",
    "INTERFACIAL",
    "HYDRATE",
    "INTER_HYDRATE",
)
# Published CHILL+ (Nguyen and Molinero, J. Phys. Chem. B 2015, 119, 9369)
# has one interfacial class and one liquid class. reCubic and reHex are not
# in that table, so they are not folded into cubic or hexagonal.
PUBLISHED_TO_OVITO = {
    0: 2,
    1: 1,
    2: 0,
    3: 3,
    4: 4,
    5: 5,
    6: 0,
}


def parse_frames(path):
    frames = []
    current = None
    for raw in open(path, encoding="utf-8"):
        line = raw.strip()
        if line == "frame":
            current = {"atoms": [], "rows": [], "labels": []}
            continue
        if current is None:
            continue
        if line == "endframe":
            frames.append(current)
            current = None
            continue
        parts = line.split()
        key = parts[0]
        if key == "step":
            current["step"] = int(parts[1])
        elif key == "cutoff":
            current["cutoff"] = float(parts[1])
        elif key == "box":
            current["box"] = [float(x) for x in parts[1:]]
        elif key == "pbc":
            current["pbc"] = [int(x) == 0 for x in parts[1:]]
        elif key == "host":
            current["host"] = int(parts[1])
        elif key == "atoms":
            current["n"] = int(parts[1])
        elif key == "atom":
            current["atoms"].append(
                (int(parts[1]), int(parts[2]), float(parts[3]), float(parts[4]), float(parts[5]))
            )
        elif key == "row":
            ids = [int(x) for x in parts[1:]]
            current["rows"].append((ids[0], set(ids[1:])))
        elif key == "label":
            current["labels"].append(
                (
                    int(parts[1]),
                    int(parts[2]),
                    float(parts[3]),
                    float(parts[4]),
                    float(parts[5]),
                    int(parts[6]),
                )
            )
    return frames


def host_sets(frame):
    return {central: neigh for central, neigh in frame["rows"]}


def freud_sets(frame):
    atoms = frame["atoms"] if frame["atoms"] else [
        (item[0], item[1], item[2], item[3], item[4]) for item in frame["labels"]
    ]
    ids = [item[0] for item in atoms]
    points = np.array([[item[2], item[3], item[4]] for item in atoms], dtype=np.float64)
    xlo, ylo, zlo, lx, ly, lz, xy, xz, yz = frame["box"]
    box = freud.box.Box(Lx=lx, Ly=ly, Lz=lz, xy=xy, xz=xz, yz=yz)
    shifted = points - np.array([xlo, ylo, zlo])
    shifted[:, 0] -= 0.5 * (lx + xy + xz)
    shifted[:, 1] -= 0.5 * (ly + yz)
    shifted[:, 2] -= 0.5 * lz
    wrapped = box.wrap(shifted)
    query = freud.locality.AABBQuery(box, wrapped)
    result = query.query(wrapped, {"r_max": frame["cutoff"], "exclude_ii": True}).toNeighborList()
    sets = {atom_id: set() for atom_id in ids}
    for query_index, point_index in zip(result.query_point_indices, result.point_indices):
        sets[ids[int(query_index)]].add(ids[int(point_index)])
    return sets


def ovito_cell(frame, data, points):
    xlo, ylo, zlo, lx, ly, lz, xy, xz, yz = frame["box"]
    cell = SimulationCell(pbc=tuple(frame["pbc"]))
    matrix = np.array(
        [[lx, xy, xz, xlo], [0.0, ly, yz, ylo], [0.0, 0.0, lz, zlo]],
        dtype=np.float64,
    )
    cell[:, :] = matrix
    data.objects.append(cell)
    particles = Particles()
    particles.create_property("Position", data=points)
    data.objects.append(particles)
    return particles


def ovito_sets(frame):
    atoms = frame["atoms"]
    ids = [item[0] for item in atoms]
    points = np.array([[item[2], item[3], item[4]] for item in atoms], dtype=np.float64)
    data = DataCollection()
    ovito_cell(frame, data, points)
    finder = CutoffNeighborFinder(frame["cutoff"], data)
    sets = {}
    for index, atom_id in enumerate(ids):
        sets[atom_id] = {ids[neigh.index] for neigh in finder.find(index)}
    return sets


def frame_points(frame):
    rows = frame["atoms"] if frame["atoms"] else [
        (item[0], item[1], item[2], item[3], item[4]) for item in frame["labels"]
    ]
    points = np.array([[item[2], item[3], item[4]] for item in rows], dtype=np.float64)
    return rows, points


def minimum_image(delta, box):
    xlo, ylo, zlo, lx, ly, lz, xy, xz, yz = box
    cell = np.array([[lx, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]], dtype=np.float64)
    frac = np.linalg.solve(cell.T, delta.reshape(-1, 3).T).T
    frac -= np.round(frac)
    return (cell.T @ frac.T).T.reshape(delta.shape)


def boundary_proof(frame):
    _rows, points = frame_points(frame)
    cutoff = frame["cutoff"]
    count = len(points)
    nearest = 1.0e300
    close = 0
    width = 256
    for start in range(0, count, width):
        block = points[start : start + width]
        delta = block[:, None, :] - points[None, :, :]
        delta = minimum_image(delta, frame["box"])
        distance = np.sqrt(np.sum(delta * delta, axis=2))
        for local, global_index in enumerate(range(start, start + len(block))):
            distance[local, global_index] = np.inf
        gap = np.abs(distance - cutoff)
        nearest = min(nearest, float(np.min(gap)))
        close += int(np.sum(gap < 1.0e-6))
    pairs = count * (count - 1) // 2
    close //= 2
    print(
        "boundary step %d atoms %d rc %.17g min_abs_r_minus_rc %.17g pairs_within_1e-6 %d pairs %d"
        % (frame["step"], count, cutoff, nearest, close, pairs),
        flush=True,
    )
    return close


def compare_sets(frame, reference, name):
    host = host_sets(frame)
    bad = 0
    for atom_id, host_set in host.items():
        if host_set != reference.get(atom_id, set()):
            bad += 1
    print(
        "step %d atoms %d %s_set_mismatches %d" % (frame["step"], len(host), name, bad),
        flush=True,
    )
    return bad


def ylm3(delta):
    radius = np.linalg.norm(delta)
    direction = delta / radius
    cosine = direction[2]
    sine = np.sqrt(max(0.0, 1.0 - cosine * cosine))
    phi = np.arctan2(direction[1], direction[0])

    def phase(order):
        return np.exp(1j * order * phi)

    return np.array(
        [
            -np.sqrt(35.0 / np.pi) / 8.0 * phase(-3) * sine ** 3,
            np.sqrt(105.0 / (2.0 * np.pi)) / 4.0 * phase(-2) * sine ** 2 * cosine,
            -np.sqrt(21.0 / np.pi) / 8.0 * phase(-1) * sine * (5.0 * cosine * cosine - 1.0),
            np.sqrt(7.0 / np.pi) / 4.0 * (5.0 * cosine ** 3 - 3.0 * cosine),
            np.sqrt(21.0 / np.pi) / 8.0 * phase(1) * sine * (5.0 * cosine * cosine - 1.0),
            np.sqrt(105.0 / (2.0 * np.pi)) / 4.0 * phase(2) * sine ** 2 * cosine,
            -np.sqrt(35.0 / np.pi) / 8.0 * phase(3) * sine ** 3,
        ]
    )


def bond_counts(points, box, index, cutoff):
    count = len(points)
    qlm = np.zeros((count, 7), dtype=complex)
    neighbours = []
    for atom in range(count):
        delta = minimum_image(points - points[atom], box)
        distance = np.sqrt(np.sum(delta * delta, axis=1))
        distance[atom] = 1.0e300
        order = np.argsort(distance)[:4]
        if distance[order[-1]] >= cutoff:
            order = np.where(distance < cutoff)[0][:4]
        neighbours.append(order)
        total = np.zeros(7, dtype=complex)
        for other in order:
            total += ylm3(delta[other])
        if len(order):
            total /= len(order)
        qlm[atom] = total
    delta = minimum_image(points - points[index], box)
    distance = np.sqrt(np.sum(delta * delta, axis=1))
    distance[index] = 1.0e300
    order = np.argsort(distance)[:4]
    eclipsed = staggered = 0
    values = []
    for other in order:
        left = qlm[index]
        right = qlm[other]
        left_norm = np.sqrt(np.vdot(left, left).real)
        right_norm = np.sqrt(np.vdot(right, right).real)
        value = (np.vdot(left, right) / (left_norm * right_norm)).real
        values.append(float(value))
        if value <= -0.8:
            staggered += 1
        elif -0.35 <= value <= 0.25:
            eclipsed += 1
    return eclipsed, staggered, [float(distance[item]) for item in order], values


def full_shells(points, box, cutoff):
    xlo, ylo, zlo, lx, ly, lz, xy, xz, yz = box
    cell = np.array([[lx, xy, xz], [0.0, ly, yz], [0.0, 0.0, lz]], dtype=np.float64)
    relative = points - np.array([xlo, ylo, zlo])
    fractional = np.linalg.solve(cell.T, relative.T).T
    shells = [None] * len(points)
    distances = [None] * len(points)
    vectors = [None] * len(points)
    width = 200
    for start in range(0, len(points), width):
        block = fractional[start : start + width]
        delta = block[:, None, :] - fractional[None, :, :]
        delta -= np.round(delta)
        cartesian = (cell.T @ delta.reshape(-1, 3).T).T.reshape(delta.shape)
        radius = np.sqrt(np.sum(cartesian * cartesian, axis=2))
        for local, row in enumerate(radius):
            index = start + local
            row[index] = np.inf
            chosen = np.where(row <= cutoff)[0]
            chosen = chosen[np.argsort(row[chosen])]
            shells[index] = chosen
            distances[index] = row[chosen]
            vectors[index] = cartesian[local, chosen]
    return shells, distances, vectors


def q_vectors(shells, vectors, full):
    count = len(shells)
    qlm = np.zeros((count, 7), dtype=complex)
    for index, chosen in enumerate(shells):
        limit = len(chosen) if full else min(4, len(chosen))
        if limit == 0:
            continue
        total = np.zeros(7, dtype=complex)
        for pos in range(limit):
            total += ylm3(vectors[index][pos])
        qlm[index] = total / limit
    return qlm


def eclipsed_staggered(index, shells, qlm):
    left = qlm[index]
    left_norm = np.sqrt(np.vdot(left, left).real)
    eclipsed = staggered = 0
    for other in shells[index][:4]:
        right = qlm[other]
        right_norm = np.sqrt(np.vdot(right, right).real)
        if left_norm == 0.0 or right_norm == 0.0:
            continue
        value = (np.vdot(left, right) / (left_norm * right_norm)).real
        if value <= -0.8:
            staggered += 1
        elif -0.35 <= value <= 0.25:
            eclipsed += 1
    return eclipsed, staggered


def chill_table(frame, source):
    rows, points = frame_points(frame)
    data = DataCollection()
    ovito_cell(frame, data, points)
    data.apply(ChillPlusModifier(cutoff=frame["cutoff"]))
    found = np.array(data.particles["Structure Type"], dtype=int)
    table = np.zeros((len(DSEAMS_NAMES), len(OVITO_NAMES)), dtype=int)
    codes = [item[5] for item in frame["labels"]]
    for code, label in zip(codes, found):
        table[code, int(label)] += 1
    print("chill_confusion step %d atoms %d source %s" % (frame["step"], len(codes), source), flush=True)
    print("category " + " ".join(OVITO_NAMES), flush=True)
    off = []
    for row, name in enumerate(DSEAMS_NAMES):
        print(name + " " + " ".join(str(int(value)) for value in table[row]), flush=True)
    for index, code in enumerate(codes):
        expect = PUBLISHED_TO_OVITO.get(code)
        if expect is None or int(found[index]) != expect:
            off.append(index)
    print("ovito_off_diagonal %d" % len(off), flush=True)
    shells, distances, vectors = full_shells(points, frame["box"], frame["cutoff"])
    near = q_vectors(shells, vectors, False)
    wide = q_vectors(shells, vectors, True)
    groups = {"A": [], "B": [], "C": []}
    ids = [item[0] for item in rows]
    for index in off:
        shell = shells[index]
        if len(shell) > 4:
            group = "A"
        else:
            near_counts = eclipsed_staggered(index, shells, near)
            wide_counts = eclipsed_staggered(index, shells, wide)
            group = "B" if near_counts != wide_counts else "C"
        groups[group].append(index)
        shell_text = " ".join(
            "%d %.17g" % (ids[int(other)], float(distances[index][pos]))
            for pos, other in enumerate(shell)
        )
        print(
            "group %s source %s atom %d code %s ovito %s n %d shell %s"
            % (
                group,
                source,
                ids[index],
                DSEAMS_NAMES[codes[index]],
                OVITO_NAMES[int(found[index])],
                len(shell),
                shell_text,
            ),
            flush=True,
        )
    print(
        "group_counts source %s A %d B %d C %d" % (source, len(groups["A"]), len(groups["B"]), len(groups["C"])),
        flush=True,
    )
    return groups


def print_repo(label, repo):
    if not repo or not Path(repo).exists():
        print("%s absent" % label, flush=True)
        return
    head = subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD"], text=True).strip()
    tree = subprocess.check_output(["git", "-C", repo, "rev-parse", "HEAD^{tree}"], text=True).strip()
    print("HEAD_%s %s" % (label, head), flush=True)
    print("TREE_%s %s" % (label, tree), flush=True)


def main():
    if freud is None:
        raise SystemExit("freud is not installed")
    print_repo("seams", os.environ.get("SRC_SEAMS", ""))
    print_repo("lammps", os.environ.get("SRC_LAMMPS", ""))
    print("freud %s" % freud.__version__, flush=True)
    print("numpy %s" % np.__version__, flush=True)
    print("python %s" % sys.version.split()[0], flush=True)
    if ovito is None:
        print("ovito unavailable", flush=True)
    else:
        print("ovito %s" % ".".join(str(part) for part in ovito.version), flush=True)
    print("rings_equivalent: none", flush=True)
    print("cages_equivalent: none", flush=True)
    print(
        "rings_note: OVITO FindRingsModifier walks a bond mesh and does not compute the primitive ring census",
        flush=True,
    )
    print(
        "chill_map: cubic->CUBIC hexagonal->HEXAGONAL interfacial->INTERFACIAL clathrate->HYDRATE interClathrate->INTER_HYDRATE water->OTHER unclassified->OTHER",
        flush=True,
    )
    print(
        "chill_map_note: Nguyen and Molinero 2015 Table 1 is not in hand, so its rows are not stated here. water and unclassified share OTHER. reCubic and reHex stay separate and are not folded into cubic or hexagonal.",
        flush=True,
    )
    print(
        "group_A_note: unresolved pending Nguyen and Molinero 2015 Table 1 and Eq. 1. More than 4 neighbours within the cutoff, so OVITO returns OTHER and d-SEAMS classifies its 4 nearest.",
        flush=True,
    )
    print(
        "group_B_note: unresolved pending Nguyen and Molinero 2015 Table 1 and Eq. 1. Exactly 4 neighbours, and one neighbour q3 uses the full shell in OVITO and the 4 nearest in d-SEAMS.",
        flush=True,
    )
    print(
        "group_C_note: conditional on the Moore 2010 interfacial-neighbour test, pending Nguyen and Molinero 2015 Table 1. d-SEAMS applies that test and OVITO ChillPlus skips it. The minimal case is mW atom 775.",
        flush=True,
    )
    host_bad = 0
    host_frames = 0
    chill_frames = 0
    totals = {"A": 0, "B": 0, "C": 0}
    for path in sys.argv[1:]:
        print("file %s" % path, flush=True)
        source = Path(path).parent.name
        for frame in parse_frames(path):
            boundary_proof(frame)
            if frame.get("host") == 1:
                host_frames += 1
                host_bad += compare_sets(frame, freud_sets(frame), "freud")
                if ovito is not None:
                    host_bad += compare_sets(frame, ovito_sets(frame), "ovito")
            if frame["labels"]:
                chill_frames += 1
                if ChillPlusModifier is None:
                    print("chill_equivalent: unavailable", flush=True)
                else:
                    groups = chill_table(frame, source)
                    for key in totals:
                        totals[key] += len(groups[key])
    print("host_frames %d host_set_mismatches %d" % (host_frames, host_bad), flush=True)
    print("chill_frames %d" % chill_frames, flush=True)
    print(
        "group_totals A %d B %d C %d" % (totals["A"], totals["B"], totals["C"]),
        flush=True,
    )
    if host_frames < 1 or chill_frames < 1:
        raise SystemExit("missing host or chill frames")
    if host_bad:
        raise SystemExit("neighbour set mismatches %d" % host_bad)
    print("exit: 0", flush=True)


if __name__ == "__main__":
    main()
