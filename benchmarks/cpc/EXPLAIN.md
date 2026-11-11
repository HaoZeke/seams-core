# Liquid labels and cage times

Both answers use the release build, one thread on core 0, and the
512-molecule mW liquid at 215 K or the ideal ice Ic ladder. Nothing
in the CHILL or CHILL+ thresholds was changed.

## Liquid agreement

Main CHILL+ at 7f310225 agrees with OVITO 3.17.2 on 439 of 512
molecules (0.857). The published rule at e8610637 agrees on 479
(0.936). The ovito-chill+ rule at that commit agrees on all 512.
The ovito-chill+ labels match OVITO's ChillPlusModifier on every
molecule, so that rule is the OVITO behaviour on this frame.
`records/explain/liquid-steps.tsv` names the step for every molecule
that still disagrees. No c_ij value lands on -0.8, -0.35 or 0.25, so
the open versus closed edge does not decide a label here.

Main and the published rule compute the same four-nearest c_ij on all
512 molecules. They differ on 40 molecules, every one of them with a
cutoff shell other than four (27 of shell 5, 12 of shell 6, one of
shell 7). Main classifies the four nearest bonds. The published rule
requires the cutoff shell itself to hold four neighbours, Nguyen and
Molinero Table 1 footnote a, and labels the other shells water. That
is the rule stated for d-SEAMS in the paper: four nearest neighbours
inside the cutoff, while OVITO scores a molecule only when the cutoff
shell has four neighbours (`paper/tex/rg_main.tex` lines 491 and 493).
Main follows the paper. The shell gate is a definition difference, not
a bug in main.

The remaining 33 molecules, all with a shell of four, are where the
published rule still differs from OVITO. Twenty-seven of them change
bond-class counts because q_3 of a neighbour is averaged over a
different set: the published rule uses the four nearest neighbours in
the shell, and OVITO averages q_3 over the whole shell
(`chill_plus_scheme.hpp` lines 28-36 at e8610637). The centre has four
neighbours, but a partner often does not, so c_ij moves. The last six
keep the same staggered and eclipsed counts. Their labels differ
because the published rule still applies the Table 1 neighbour test
and ovito-chill+ does not: two staggered bonds, or three staggered and
no eclipsed bond, are interfacial without looking at the neighbour.
Those six are the same definition, not a bug.

Verdict: no BUG. No product change.

## Cage times

The paper's double-diamond median of 2.61 ms, whole path of 5.28 ms,
and exponents 0.761 and 0.889 are the log-log least-squares fit in
`paper/tex/figures/stats.py` to `dseams2-ic-*.txt`. Those logs cover
216, 512, 1000, 1728, 2744 and 4096 molecules. The 1000-molecule log
is thirty calls after five warm-ups. The machine line recorded with
them is an Intel Xeon, model 207, and GCC 15.2.0. The tag those
sources belong to is v2.12.0, commit 0fc9ec49.

On this machine the same protocol as the rest of this study, one
warm-up and six repeats, rebuilds v2.12.0 and returns 2.591 ms for
the double-diamond stage and 5.120 ms for the whole path at 1000
molecules. The exponents on the six paper sizes are 0.758 and 0.884.
The current sources, which match commit 7e9a2f4c, return 3.270 ms and
5.385 ms, with exponents 0.823 and 0.900. Leaving 2744 out of the
paper fit only moves 0.761 to 0.764, so the fit range is not the
shift. The hardware rerun of v2.12.0 reproduces the paper's
double-diamond time. The slower current stage is the cage code after
that tag: `4b62df96` packs the neighbour rows for one ring search,
`df2a2be3` keeps each atom's rings in one list, and `7e9a2f4c` stores
neighbour bonds in sorted rows. The rings at 1000 molecules got
faster and the double-diamond search got slower, and the whole path
moved only from 5.12 ms to 5.39 ms.

A release build of tag v2.12.0 on this machine, one thread and the
same six ice Ic sizes, recovers the paper's double-diamond time of
2.61 ms and its exponents. The current tree is slower in that stage,
3.27 ms, because the ring and cage packing changed after the tag;
the whole path stays near 5.3 ms because the ring search got faster
by about as much as the double-diamond search got slower. The
exponent shift has the same cause. It is not the processor, the
compiler, or the choice of sizes in the log-log fit.
