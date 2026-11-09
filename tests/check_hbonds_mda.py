# Copyright (c) d-SEAMS developers
"""Compare seams hbonds with MDAnalysis on one frame.

usage: check_hbonds_mda.py SEAMS LAMMPSTRJ

d-SEAMS prints the number of O-O pairs that pass the test.
MDAnalysis HydrogenBondAnalysis prints donor-H-acceptor rows.
One O-O pair can be two of those rows when each molecule donates.
"""

import subprocess
import sys

import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import (
    HydrogenBondAnalysis,
)

seams, path = sys.argv[1], sys.argv[2]
proc = subprocess.run(
    [seams, "hbonds", path, "--frame", "1", "--last", "1", "--type", "2",
     "--htype", "1", "--cutoff", "3.5", "--hdist", "2.5", "--hangle", "30"],
    capture_output=True, text=True)
sys.stdout.write(proc.stdout)
if proc.returncode != 0:
    sys.stderr.write(proc.stderr)
    sys.exit(proc.returncode)

u = mda.Universe(path, format="LAMMPSDUMP")
# LAMMPS type 2 is oxygen in these frames, type 1 is hydrogen.
hb = HydrogenBondAnalysis(
    u,
    donors_sel="type 2",
    hydrogens_sel="type 1",
    acceptors_sel="type 2",
    d_a_cutoff=3.0,
    d_h_a_angle_cutoff=150.0,
)
hb.run()
rows = hb.results.hbonds
print("MDAnalysis directed bonds", len(rows))
for row in rows:
    donor = int(row[1])
    hydrogen = int(row[2])
    acceptor = int(row[3])
    dist = float(row[4])
    angle = float(row[5])
    print("D", donor, "H", hydrogen, "A", acceptor,
          "dist", round(dist, 3), "angle", round(angle, 2))
