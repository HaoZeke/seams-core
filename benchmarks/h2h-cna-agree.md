# Adaptive CNA against CHILL+

The three correctness rows are ice Ih, ice Ic, and the same Ih cell after a 1.5 A Gaussian shift (PCG64 seed 20261113), each 4096 molecules. CHILL+ is d-SEAMS ovito mode. Adaptive CNA is OVITO CommonNeighborAnalysisModifier mode AdaptiveCutoff, the call in scripts/h2h_summary_extra.py.

A periodic FCC copper cell, lattice constant 3.615 A, 4x4x4 cubic cells, 256 atoms, goes through that same modifier: 256 atoms are FCC. The ice result is not a pipeline that labels everything Other.

Adaptive CNA writes Structure Type Other, FCC, HCP, BCC, or ICO (ovito/modifiers/__init__.pyi). Those signatures are close-packed or bcc shells. CHILL+ writes a water class from the staggered and eclipsed counts of a four-neighbor shell (src/include/internal/bop.hpp). A four-neighbor shell matches none of the four crystal types, so the structure type is Other.

| Phase | CHILL+ | Adaptive CNA | Coordination | Atoms | Cause |
| --- | --- | --- | --- | --- | --- |
| Ih | hexagonal_ice | Other | 4 | 4096 | CHILL+ row on a 4-neighbor shell; adaptive CNA has no tetrahedral type |
| Ic | cubic_ice | Other | 4 | 4096 | CHILL+ row on a 4-neighbor shell; adaptive CNA has no tetrahedral type |
| liquid | hexagonal_ice | Other | 4 | 1 | CHILL+ row on a 4-neighbor shell; adaptive CNA has no tetrahedral type |
| liquid | hydrate | Other | 4 | 49 | CHILL+ row on a 4-neighbor shell; adaptive CNA has no tetrahedral type |
| liquid | interfacial_hydrate | Other | 4 | 168 | CHILL+ row on a 4-neighbor shell; adaptive CNA has no tetrahedral type |
| liquid | interfacial_ice | Other | 4 | 4 | CHILL+ row on a 4-neighbor shell; adaptive CNA has no tetrahedral type |
| liquid | other | Other | 0 | 8 | both unstructured, or the boolean labels match |
| liquid | other | Other | 1 | 63 | both unstructured, or the boolean labels match |
| liquid | other | Other | 2 | 261 | both unstructured, or the boolean labels match |
| liquid | other | Other | 3 | 566 | both unstructured, or the boolean labels match |
| liquid | other | Other | 4 | 634 | both unstructured, or the boolean labels match |
| liquid | other | Other | 5 | 830 | both unstructured, or the boolean labels match |
| liquid | other | Other | 6 | 723 | both unstructured, or the boolean labels match |
| liquid | other | Other | 7 | 421 | both unstructured, or the boolean labels match |
| liquid | other | Other | 8 | 221 | both unstructured, or the boolean labels match |
| liquid | other | Other | 9 | 108 | both unstructured, or the boolean labels match |
| liquid | other | Other | 10 | 29 | both unstructured, or the boolean labels match |
| liquid | other | Other | 11 | 5 | both unstructured, or the boolean labels match |
| liquid | other | Other | 12 | 4 | both unstructured, or the boolean labels match |
| liquid | other | Other | 13 | 1 | both unstructured, or the boolean labels match |

On Ih every atom is hexagonal ice, coordination 4, adaptive CNA Other. That is 3 staggered bonds and 1 eclipsed bond. On Ic every atom is cubic ice, coordination 4, adaptive CNA Other. That is 4 staggered bonds. The 222 liquid disagreements are the waters that still have exactly four neighbors inside 3.5 A and still match one CHILL+ row after the shift: 1 hexagonal, 49 hydrate, 168 interfacial hydrate, and 4 interfacial ice. Every one of them is adaptive CNA Other. Another 634 waters also have four neighbors and stay CHILL+ other, because the staggered and eclipsed counts match no row; those agree with adaptive CNA Other. Four liquid atoms have 12 neighbors inside 3.5 A and are still Other. Coordination by itself is not an FCC signature. The FCC control is a copper lattice; the shifted ice is not.

Per-atom rows are benchmarks/h2h-cna-atoms.csv.
