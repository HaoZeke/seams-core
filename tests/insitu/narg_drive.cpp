// Copyright (c) d-SEAMS developers
// The LAMMPS compute base rejects a short command before the body runs.
// This driver calls the same guard with two arguments so that arm is taken.

#include "compute_dseams.h"
#include "exceptions.h"
#include "lammps.h"

#include <mpi.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  LAMMPS_NS::LAMMPS *lmp = new LAMMPS_NS::LAMMPS(0, nullptr, MPI_COMM_WORLD);
  int fired = 0;
  try {
    LAMMPS_NS::require_dseams_narg(lmp, 2);
  } catch (LAMMPS_NS::LAMMPSException &) {
    fired = 1;
  }
  LAMMPS_NS::require_dseams_narg(lmp, 4);
  delete lmp;
  MPI_Finalize();
  return fired == 1 ? 0 : 1;
}
