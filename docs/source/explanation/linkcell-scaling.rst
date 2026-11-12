Linkcell strong scaling
=======================

An orthorhombic frame is cut by linkcell's weighted recursive coordinate
bisection of the periodic cell grid. The balance bound is the tightest a
cell-aligned cut can promise: one plus the part count times the heaviest
cell, over the total weight. The first call of an analysis weights cells
by how many points they hold. A later call of that analysis weights them
by the work measured on the first call.

The halo is painted in the same walk that visits the cells a rank owns.
A point in an owned cell stays owned. A point in a painted cell that
another rank owns is a ghost. The cut, the cell of each point, and the
neighbour rows of an unchanged frame are kept and refilled in place.

A skewed cell does not use that mesh. It keeps the slab, pencil, or
block cut.

Session
-------

One sitting, 12 November 2026. ``nproc`` was read with ``OMP_NUM_THREADS``
unset and reported 8 cores. The compiler is g++ 15.2.0. The library is a
Meson release build. Eight samples follow one warm-up. The median of an
even count is the mean of the two central samples. Speedup is against
one rank and one thread. The four operations are the neighbour rows,
CHILL+, the Steinhardt q6, and primitive rings of depth 6. The ice is
eight molecules in a cubic cell. No host name is recorded.

Eight threads share one cell mesh, so the neighbour walk does not reach
an eightfold speedup. Eight ranks on the larger frame do: q6 is 8.22
times the one-rank time, and the ring census is 7.31 times.

OpenMP, 110592 molecules
------------------------

======= ======== ======= ======== ======= ======= ========= =======
threads neigh ms speedup chill ms speedup q6 ms   speedup   rings ms
======= ======== ======= ======== ======= ======= ========= =======
1       15.4     1.00    33.2     1.00    39.8    1.00      75.1
2       9.6      1.61    18.9     1.76    22.2    1.79      46.6
4       9.6      1.61    13.9     2.39    14.7    2.71      32.6
8       8.3      1.87    11.1     2.99    11.0    3.62      23.4
======= ======== ======= ======== ======= ======= ========= =======

The rings column above is the median in milliseconds. Its speedups are
1.00, 1.61, 2.31 and 3.21.

OpenMP, 884736 molecules
------------------------

======= ======== ======= ======== ======= ======= ========= =======
threads neigh ms speedup chill ms speedup q6 ms   speedup   rings ms
======= ======== ======= ======== ======= ======= ========= =======
1       166.4    1.00    340.2    1.00    556.7   1.00      618.4
2       113.9    1.46    217.9    1.56    361.9   1.54      397.2
4       94.2     1.77    149.8    2.27    278.5   2.00      284.3
8       94.7     1.76    133.6    2.55    250.8   2.22      223.4
======= ======== ======= ======== ======= ======= ========= =======

Ring speedups are 1.00, 1.56, 2.17 and 2.77.

MPI, one thread per rank
------------------------

========= ===== ======== ======= ======== ======= ===== ======= ======== =======
molecules ranks neigh ms speedup chill ms speedup q6 ms speedup rings ms speedup
========= ===== ======== ======= ======== ======= ===== ======= ======== =======
110592    1     15.4     1.00    33.2     1.00    39.8  1.00    75.1     1.00
110592    2     18.7     0.83    29.6     1.12    32.9  1.21    39.9     1.88
110592    4     10.3     1.49    17.3     1.92    24.0  1.66    21.4     3.51
110592    8     4.1      3.76    7.6      4.38    8.1   4.91    10.7     7.02
884736    1     166.4    1.00    340.2    1.00    556.7 1.00    618.4    1.00
884736    2     127.8    1.30    238.0    1.43    328.2 1.70    336.3    1.84
884736    4     62.2     2.68    115.9    2.94    156.7 3.55    168.0    3.68
884736    8     26.4     6.30    59.3     5.74    67.8  8.22    84.6     7.31
========= ===== ======== ======= ======== ======= ===== ======= ======== =======
