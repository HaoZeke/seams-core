# Liquid labels and cage times

Both answers use the release build, one thread on core 0. The liquid
is the 512-molecule mW frame at 215 K. The cage ladder is ideal ice Ic.
No CHILL or CHILL+ threshold was retuned.

## Liquid agreement

Main CHILL+ at 7f310225 agrees with OVITO 3.17.2 on 439 of 512
molecules (0.857). The published rule at e8610637 agrees on 479
(0.936). `ovito-chill+` at that commit agrees on all 512, and those
labels match OVITO's ChillPlusModifier on this frame.
`records/explain/liquid-steps.tsv` names every molecule that disagrees.

Forty of the 73 main disagreements are a bug. `assignIceTypePlus` in
`src/bop.cpp` reads `c_ij.size()` (line 1118) and accepts a crystal row
only when that count is 4 (line 1135). `classifyBondsLeftFp64` has
already kept at most four neighbours, the CHILL+ coordination, before
it writes `c_ij`. A molecule with five, six, or seven neighbours inside
3.5 angstrom still has four recorded bonds, so the gate lets it through.
On this liquid that is 27 molecules of shell 5, 12 of shell 6, and one
of shell 7. Main labels them hydrate, interfacial hydrate, or
interfacial. OVITO labels every one of them other, and OVITO labels
other for every shell that is not exactly four. The comment above the
gate says the bound is the recorded bonds, not the neighbour list,
which is how the extra shell neighbours disappear. Table 1 of Nguyen
and Molinero requires the cutoff shell itself to contain four
neighbours, and the neighbour column in `src/include/internal/bop.hpp`
lists 4 for every crystal phase. The CHILL+ pull-request series,
including the published rule at e8610637, applies that shell count.
The agreement of that fixed code on this liquid is pending the CHILL+
pr2 tip. That number is not invented here.

The other 33 disagreements all have a shell of exactly four. They are
the same 33 molecules on which the published rule still differs from
OVITO, 479 of 512. Equation 1 averages q_3 over the four nearest
neighbours. OVITO 3.17.2 averages q_3 over the whole cutoff shell.
That is a definition difference, not a second copy of the gate bug.
No c_ij on this frame sits on -0.8, -0.35, or 0.25.

## Cage times

The paper still quotes the v2.12.0 measurement: double-diamond 2.61 ms
and exponent 0.761 on ice Ic, from the log-log fit in
`paper/tex/figures/stats.py` over 216, 512, 1000, 1728, 2744, and 4096
molecules. Rebuilding tag v2.12.0 (`0fc9ec49`) on this machine, one
warm-up and six repeats, returns 2.591 ms and exponent 0.758.

The current sources return 3.270 ms at 1000 molecules and a
double-diamond exponent of 0.823 on those six sizes (0.826 on the
five sizes that omit 2744). That slower stage is a regression
introduced in `df2a2be3` (`perf(cages): keep the rings of each atom in
one list`), which rewrites cage packing between 2.12.0 and 2.13.0.
Rings got faster and the double-diamond search got slower, so the
whole path stays near 5.39 ms. The machine is the same class the paper
records, an Intel Xeon with GCC 15.2.0, so the shift is that commit,
not the processor.

A C++ fix is on `perf/ddc-regression-fix` at `2c0db67f`. It is under
independent verify. The reported figures are about 2.74 ms at 1000
molecules and exponent 0.799. Those figures are recorded here as
pending that verify. They were not rerun in this tree.
