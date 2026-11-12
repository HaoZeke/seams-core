Cutoff neighbours, local $q_4$ and $q_6$, the diamond
common-neighbour labels, metal CNA, primitive rings and local ice
symmetry, against mdapy, pyscal3, vesin, IRA and SQQ.

The statistics are `scripts/h2h_protocol.py`: one warm-up, eight timed
calls in one process, the median of the even sample, and type-7
quartiles. One thread is pinned to cpu 0. Eight threads are pinned to
cpus 0-7. The allowed cpu list is recorded. The host name is not.
This session saw 8 cpus, `0-7`, on an Intel Xeon. d-SEAMS is 2.14.0.
The symmetry backend in this build is the in-tree Kabsch kernel;
libira is not linked.

| code | version | what was timed |
| --- | --- | --- |
| mdapy | 2.0.0 | neighbours, q4, q6, diamond, metal CNA, rings |
| pyscal3 | 4.1.0 | neighbours, q4, q6, diamond, metal CNA |
| vesin | 0.6.2 | neighbour lists |
| IRA | 2.2.0 | shape match of the 17-atom ice shell |
| SQQ | 0.5.8 | shortest-path chordless rings of size 4 to 6 |


IRA is the conda-forge package `ira` (`import ira_mod`). The PyPI
package also named `ira` is version 0.0 and has no `ira_mod`, so it is
not the shape-matching library. SQQ 0.5.8 installs from PyPI.

Frames come from `scripts/gen_ideal_ice.py` and are not stored. The
oxygen bond is 2.75 Å and the cutoff is 3.5 Å. Ice Ic is the 8-atom
cubic cell. Ice Ih is the 8-atom orthorhombic cell made from two
hexagonal primitives, so the box has no tilt. The molecule counts are
128, 4096 and 110592. Two frames already in the tree are included.
`input/traj/dump-l.lammpstrj` frame 1 is 512 atoms of type 1. Its
number density is 0.311 Å^{-3} and the median number of neighbours
inside 3.5 Å is 56, so it is not liquid water; the domain test calls
it the liquid frame. `input/traj/exampleTraj.lammpstrj` contributes
its 250 type-2 oxygens, the mixed TIP4P example.

Agreement is per atom. Neighbour rows are the same sets. $q_l$ agrees
when every atom is within $10^{-6}$. Diamond labels are compared after
mapping d-SEAMS cubic (1) and hexagonal (0) onto the mdapy and pyscal3
codes 1 and 4. Metal CNA is a different label set and is reported
against the diamond label so the mismatch is visible. Rings agree when
the size census matches. Local symmetry agrees when the IRA template
with the smaller Hausdorff distance is the same ice label, and a match
worse than 0.2 Å is left unclassified.

`snakemake literature_h2h` regenerates `benchmarks/h2h-literature.csv`
and the raw text log `benchmarks/h2h-literature.txt`. Set `H2H_PYTHON`
to an interpreter with the five packages, and `SEAMS_CORE_LIB` to the
release `libyodaLib.so`. `benchmarks/h2h_literature.py` is the one-call
Airspeed Velocity arm of the same commands, one class for each code,
at one thread and at eight.

**Neighbours.** The cutoff list matches vesin, mdapy and pyscal3 on
every frame (agreement 1). On ice Ic of 110592 molecules the medians
are 31.603 ms for d-SEAMS, 51.939 ms for vesin, 45.102 ms for mdapy and
225.169 ms for pyscal3, at one thread. At eight threads the same list
is 10.358 ms against vesin at 52.525 ms. Ice Ih of the same size is
25.773 ms against vesin at 47.958 ms. The 128-molecule cells and the
512-atom dense frame used to miss the one-cell stencil, because that
stencil required four cells on every edge. Two cells are enough: the
other image is at least half an edge away, and half an edge is at
least the cutoff, so the strict test keeps a single image. After that
change the 128-molecule ice Ic list is 0.059 ms against vesin at
0.095 ms, and the dense frame is 1.384 ms against mdapy at 1.635 ms
and vesin at 1.992 ms, still with the same neighbour sets.

**$q_4$ and $q_6$.** Local Steinhardt parameters on the cutoff list.
On ideal ice the first shell is a perfect tetrahedron, so Ih and Ic
share the same $q_l$. Every compared value is within $10^{-6}$; the
largest absolute gap on the 110592-molecule cells is a few times
$10^{-16}$. At one thread on ice Ic of that size, $q_4$ is 57.042 ms
against mdapy 62.823 ms and pyscal3 90.901 ms, and $q_6$ is 48.591 ms
against mdapy 75.009 ms and pyscal3 102.943 ms. At eight threads
$q_6$ is 13.883 ms against mdapy 95.834 ms.

**Diamond and metal CNA.** Diamond labels agree on the ideal cells:
every ice Ic molecule is cubic, every ice Ih molecule is hexagonal.
They also agree on the dense frame and on the TIP4P oxygens: every
atom is other. Metal CNA (fcc, hcp, bcc, icosahedral) labels every ice
molecule other. That is the correct CNA answer. Those signatures are
not ice classes. The diamond common-neighbour label is the one that
names cubic and hexagonal ice, and it matches. On the dense frame at
one thread, mdapy's diamond identifier is faster (0.576 ms against
1.632 ms) and returns the same all-other labels. Its time is below the
cost of the 3.5 Å list, so it is not walking that list; d-SEAMS builds
the list and then rejects every atom at the tetrahedron test. At eight
threads the same call is 1.470 ms against mdapy at 10.356 ms.

**Rings.** Primitive rings of size at most 6. On ice Ic, and on ice Ih
at 4096 and 110592 molecules, the census is two six-rings per molecule
and no smaller ring. d-SEAMS, mdapy and SQQ agree. At 110592 molecules
of ice Ic and one thread the medians are 116.834 ms, mdapy 2651.418 ms
and SQQ 18672.033 ms. At eight threads d-SEAMS is 34.116 ms and mdapy
is 2670.661 ms. SQQ does not search three-rings, and one frame does
not split across its workers, so eight threads do not change its time.
The 128-molecule ice Ih cell is the exception. d-SEAMS and SQQ both
report 32 four-rings and 384 six-rings. mdapy reports 0 and 256.
Those 32 four-rings are chordless cycles of the bond graph: the
opposite vertices are graph distance 2 and there is no chord. They are
primitive. They close through the periodic boundary of the 2 by 2 by 4
orthorhombic cell. At 4096 molecules and above the extra cycles are
gone and the three codes agree. The d-SEAMS census is the more
complete primitive count on the small cell. The dense frame has median
coordination 56, and the primitive search returns 125065 three-rings,
118033 four-rings, 2986554 five-rings and 945837 six-rings (12998 ms
at one thread). mdapy and SQQ were not asked to repeat that search.

**IRA.** Each atom's 16 nearest neighbours inside 5 Å, plus the atom
itself, are matched to the ice Ih and ice Ic templates. The template
with the smaller Hausdorff distance is the label, unless that distance
exceeds 0.2 Å. On the ideal cells the chosen distance is 0 and every
label matches `seams_symmetry_classify`. At 110592 molecules of ice Ic
the classify median is 25.236 ms at one thread and 6.562 ms at eight.
IRA is 26972 ms and 27687 ms. The eight-core IRA time does not fall,
because the matcher is one permutation search per shell. Both answers
are the same ice label. The d-SEAMS number is the bond-graph kernel,
not a call into libira. IRA remains the right tool for a cluster that
is not one of these ice shells.

The raw lines, including every sample, are
`benchmarks/h2h-literature.txt`.

