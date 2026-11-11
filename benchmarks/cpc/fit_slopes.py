# Copyright (c) d-SEAMS developers
"""Least-squares slopes from the raw cage repeats.

Each series takes the type-7 median of ddc_ms at each size, then the
least-squares slope of ln(median) on ln(molecules). One stdout line
names the series, the sizes, the slope, and the input files.
"""

from __future__ import annotations

import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RECORDS = ROOT / "benchmarks" / "cpc" / "records"
OUT = RECORDS / "explain" / "slopes.txt"


def type7(values: list[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    count = len(ordered)
    if count == 0:
        raise ValueError("empty sample")
    if count == 1:
        return ordered[0]
    position = (count - 1) * probability + 1.0
    index = int(math.floor(position))
    if index >= count:
        return ordered[-1]
    if index < 1:
        return ordered[0]
    weight = position - index
    return ordered[index - 1] + weight * (ordered[index] - ordered[index - 1])


def repeats(path: Path, column: str) -> tuple[int, list[float]]:
    molecules = None
    values: list[float] = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0] == "count" and parts[1] == "atoms":
            molecules = int(parts[2])
        if len(parts) == 4 and parts[0] == "repeat" and parts[2] == column:
            values.append(float(parts[3]))
    if molecules is None or not values:
        raise ValueError(f"no {column} repeats in {path}")
    return molecules, values


def slope(points: list[tuple[int, float]]) -> float:
    xs = [math.log(n) for n, _ in points]
    ys = [math.log(t) for _, t in points]
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    sxx = sum((x - mean_x) ** 2 for x in xs)
    return sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / sxx


def rel(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def current_file(size: int) -> Path:
    if size == 2744:
        return RECORDS / "explain" / "head-ic-2744.csv"
    return RECORDS / "b-cages" / f"ic-{size}.csv"


def v212_file(size: int) -> Path:
    return RECORDS / "explain" / f"v212-ic-{size}.csv"


SERIES = (
    ("current-ddc", (216, 512, 1000, 1728, 2744, 4096), current_file),
    ("current-ddc-omit-2744", (216, 512, 1000, 1728, 4096), current_file),
    ("v212-ddc", (216, 512, 1000, 1728, 2744, 4096), v212_file),
)


def line_for(name: str, sizes: tuple[int, ...], finder) -> str:
    points = []
    files = []
    for size in sizes:
        path = finder(size)
        molecules, values = repeats(path, "ddc_ms")
        if molecules != size:
            raise ValueError(f"{path} has {molecules} molecules, expected {size}")
        points.append((molecules, type7(values, 0.5)))
        files.append(rel(path))
    fitted = slope(points)
    size_text = ",".join(str(size) for size in sizes)
    file_text = ",".join(files)
    return (
        f"series {name} sizes {size_text} slope {format(fitted, '.16g')} "
        f"files {file_text}"
    )


def main() -> None:
    lines = [line_for(name, sizes, finder) for name, sizes, finder in SERIES]
    text = "\n".join(lines) + "\n"
    OUT.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
