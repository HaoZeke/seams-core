# Liquid labels and cage times

We measure both answers on the release build, with one thread on core 0.
The liquid is the 512-molecule mW frame at 215 K. The cage ladder is
ideal ice Ic. We left every CHILL and CHILL+ threshold at its published
value.

## Liquid agreement

Main CHILL+ at 7f310225 agrees with OVITO 3.17.2 on 439 of 512
molecules (0.857). The published rule at e8610637 agrees on 479
(0.936). The ovito-chill+ option at that commit agrees on all 512.
Those labels match OVITO's ChillPlusModifier on this frame.
`records/explain/liquid-steps.tsv` names every molecule that disagrees.

Forty of the 73 main disagreements come from a size check that reads
the wrong count. `assignIceTypePlus` in `src/bop.cpp` reads
`c_ij.size()` (`src/bop.cpp:1118`). It accepts a crystal row only when
that count is 4 (`src/bop.cpp:1135`). `classifyBondsLeftFp64` has already kept at
most four neighbours, the CHILL+ coordination, before it writes
`c_ij`. A molecule with five to seven neighbours inside 3.5 angstrom
still has four recorded bonds, so the check lets it through. On this
liquid the extra shells are shell 5 for 27 molecules, shell 6 for 12
molecules, and shell 7 for one molecule. Main labels those molecules
hydrate or interfacial. OVITO labels each of them other. OVITO also
labels other whenever the shell is not exactly four. The comment above
the check says the bound is the recorded bond count. Table 1 of Nguyen
and Molinero requires the cutoff shell itself to contain four
neighbours. The neighbour column in `src/include/internal/bop.hpp`
lists 4 for every crystal phase. The CHILL+ pull-request series
applies that shell count, and the published rule at e8610637 does too.
The agreement of that corrected code on this liquid waits on the
CHILL+ pr2 tip. We are leaving that number unmeasured.

The other 33 disagreements all have a shell of exactly four. They are
the same 33 molecules on which the published rule still differs from
OVITO, which is 479 of 512. Equation 1 averages q3 over the four
nearest neighbours. OVITO 3.17.2 averages q3 over the whole cutoff
shell. That difference is the published q3 definition. No c_ij on this
frame sits on -0.8, -0.35, or 0.25.

## Cage times

The paper quotes the v2.12.0 measurement. Double-diamond time there is
2.61 ms, and the exponent is 0.761 on ice Ic. The fit is the log-log
slope in `paper/tex/figures/stats.py`. The sizes are 216, 512, 1000,
1728, 2744, 4096 molecules. We rebuilt tag v2.12.0 (`0fc9ec49`) on
this machine with one warm-up and six repeats. The type-7 median is
2.590806 ms, and the exponent is 0.758.

The current sources return a type-7 median of 3.270127 ms at 1000
molecules. The double-diamond exponent on the six sizes is 0.823.
Omitting 2744 moves that exponent to 0.826. Both builds report 2000
equators and 500 double-diamond cages. The cage counts are identical.
The slower stage follows the packing rewrite in `df2a2be3`, which keeps
the rings of each atom in one list. Ring search time fell, and the
double-diamond search time rose. The whole path stays near 5.39 ms.
The machine is the same class the paper records, an Intel Xeon with
GCC 15.2.0.

A later change on `perf/ddc-regression-fix` at `2c0db67f` is under
separate verify. The reported figures are about 2.74 ms at 1000
molecules and exponent 0.799. We did not rerun those figures in this
tree.
