# Orthorhombic sH supercell

SQQ 0.5.8 rejects the tilted genice sH dump. The cell has `xy = -lx/4` and no `xz` or `yz` tilt. Four images along `b` fill the rectangle `lx` by `4*ly` by `lz` (24.842423 by 43.028338 by 20.131291 A) with 544 molecules and no tilt. SQQ accepts that dump.

Cage counts on that cell match. The seams command and SQQ both report 48 cages of `5:12`, 16 of `5:12,6:8`, and 32 of `4:3,5:6,6:3`. That is four times the tilted unit cell (12, 4 and 8).

Eight interleaved samples after one warm-up, commit `ef38adb99d5f8236b9d6dcae0d701421fcafed99` on every line. SQQ's python cage walk does not speed up from one thread to eight (median 1190.1 ms and 1191.8 ms). The seams command is 16.1 ms and 11.9 ms. The raw log is `benchmarks/h2h-sqq-sh-ortho.txt`. The tilted dump is still rejected by `tests/test_sqq_sh.py`.
