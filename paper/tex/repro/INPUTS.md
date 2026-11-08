The campaign resolves its paths from the repository root.
SEAMS_TIP defaults to bbdir-release.
SEAMS_BASE defaults to bbdir-base.
That directory is a build of commit 1c8efa38.
docs/orgmode/paper/versions.org names that commit.
It is the d-SEAMS 1 baseline.
SEAMS_1E6 defaults to input/traj/ice-1e6.lammpstrj.
The paper-ice-1e6 task produces that frame.
The paper records 995328 molecules for it.
The cubic cell input/traj/mW_cubic.lammpstrj has 4096 atoms.
A 9 by 9 by 3 replica of that cell has 995328 atoms.
The task paper-ice-1e6 writes the frame.
A missing base directory stops the run and names the path.
A missing frame prints its own path and the task that writes it.
