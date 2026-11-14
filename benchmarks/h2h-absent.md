# Absent summary rows

The summary had 14 absent rows: pyscal3 CHILL+ at both widths, freud rings at three sizes and both widths, and ClavisCage CHILL+ at three sizes and both widths.

## pyscal3 CHILL+

pyscal3 4.1.0 __all__ has no CHILL or CHILL+ name. common_neighbor_analysis is adaptive CNA and returns fcc, hcp, bcc, ico, and others (pyscal3/descriptors.py). diamond_structure is extended CNA. Neither is Nguyen and Molinero CHILL+. PyPI's latest pyscal3 is this same 4.1.0.

The two rows stay absent. There is no equivalent call to time.

## freud rings

freud 3.6.1 modules are box, cluster, data, density, diffraction, environment, errors, interface, locality, msd, order, parallel, plot, pmft, util. There is no ring or cycle module. order.py is Steinhardt, Hexatic, SolidLiquid, Cubatic, Nematic, RotationalAutocorrelation, and ContinuousCoordination. Ramasubramani and others, Comput. Phys. Commun. 254, 107275 (2020).

Neighbor lists are already a separate timed row. A ring search written on top of that list would be a new code, not freud. The six rows stay absent.

## ClavisCage CHILL+

ClavisCage 0.3.18 imported as clavis_engine. CHILL+ is calculate_chill_cell in src/chill.rs of Zenodo 10.5281/zenodo.24112918. The six summary rows are that call, pinned and interleaved with d-SEAMS.

PyPI's package named clavis is unrelated. Its module docstring is a credential client.
