// Copyright (c) d-SEAMS developers
// SPDX-License-Identifier: MIT
#include "compute_dseams.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"
#include "utils.h"

#include "seams.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <vector>

using namespace LAMMPS_NS;

ComputeDseams::ComputeDseams(LAMMPS *lmp, int narg, char **arg)
    : Compute(lmp, narg, arg), nmax_(0), label_(nullptr), cutoff_(3.5),
      graph_(SEAMS_GRAPH_CUTOFF), rule_(SEAMS_CHILL_PUBLISHED), calls_(0), seconds_(0.0) {
  if (narg < 3) error->all(FLERR, "Illegal compute dseams command");
  for (int iarg = 3; iarg < narg; iarg += 2) {
    if (iarg + 1 >= narg) error->all(FLERR, "Illegal compute dseams command");
    if (strcmp(arg[iarg], "cutoff") == 0) {
      cutoff_ = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (!(cutoff_ > 0.0)) error->all(FLERR, "compute dseams: cutoff must be positive");
    } else if (strcmp(arg[iarg], "graph") == 0) {
      if (strcmp(arg[iarg + 1], "mutual") == 0) {
        graph_ = SEAMS_GRAPH_MUTUAL;
      } else if (strcmp(arg[iarg + 1], "cutoff") == 0) {
        graph_ = SEAMS_GRAPH_CUTOFF;
      } else {
        error->all(FLERR, "compute dseams: graph is cutoff or mutual");
      }
    } else if (strcmp(arg[iarg], "rule") == 0) {
      if (strcmp(arg[iarg + 1], "ovito") == 0) {
        rule_ = SEAMS_CHILL_OVITO;
      } else if (strcmp(arg[iarg + 1], "published") == 0) {
        rule_ = SEAMS_CHILL_PUBLISHED;
      } else {
        error->all(FLERR, "compute dseams: rule is published or ovito");
      }
    } else {
      error->all(FLERR, "Illegal compute dseams command");
    }
  }
  peratom_flag = 1;
  size_peratom_cols = 0;
}

ComputeDseams::~ComputeDseams() {
  memory->destroy(label_);
  if (comm->me == 0 && calls_ > 0) {
    char line[160];
    std::snprintf(line, sizeof(line),
                  "compute dseams calls %ld seconds %.9e per_call %.9e\n", calls_, seconds_,
                  seconds_ / static_cast<double>(calls_));
    utils::logmesg(lmp, line);
  }
}

void ComputeDseams::init() {}

void ComputeDseams::compute_peratom() {
  invoked_peratom = update->ntimestep;
  const double t0 = MPI_Wtime();
  if (atom->nmax > nmax_) {
    memory->destroy(label_);
    nmax_ = atom->nmax;
    memory->create(label_, nmax_, "dseams:label");
    vector_atom = label_;
  }
  const int nlocal = atom->nlocal;
  std::vector<tagint> tags;
  std::vector<double> xyz;
  for (int i = 0; i < nlocal; ++i) {
    label_[i] = -1.0;
    if (!(atom->mask[i] & groupbit)) continue;
    tags.push_back(atom->tag[i]);
    xyz.insert(xyz.end(), {atom->x[i][0], atom->x[i][1], atom->x[i][2]});
  }
  const int nprocs = comm->nprocs;
  const int me = comm->me;
  int mine = static_cast<int>(tags.size());
  std::vector<int> counts(nprocs), displs(nprocs), counts3(nprocs), displs3(nprocs);
  MPI_Gather(&mine, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, world);
  int total = 0;
  if (me == 0) {
    for (int p = 0; p < nprocs; ++p) {
      displs[p] = total;
      counts3[p] = 3 * counts[p];
      displs3[p] = 3 * total;
      total += counts[p];
    }
  }
  std::vector<tagint> all_tags(me == 0 ? total : 0);
  std::vector<double> all_xyz(me == 0 ? 3 * static_cast<std::size_t>(total) : 0);
  MPI_Gatherv(tags.data(), mine, MPI_LMP_TAGINT, all_tags.data(), counts.data(), displs.data(),
              MPI_LMP_TAGINT, 0, world);
  MPI_Gatherv(xyz.data(), 3 * mine, MPI_DOUBLE, all_xyz.data(), counts3.data(), displs3.data(),
              MPI_DOUBLE, 0, world);
  MPI_Bcast(&total, 1, MPI_INT, 0, world);

  // Rank 0 classifies the group in atom-ID order and broadcasts the labels.
  std::vector<tagint> sorted_tags(total);
  std::vector<int> labels(total);
  int status = SEAMS_OK;
  if (me == 0 && total > 0) {
    std::vector<int> order(total);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return all_tags[a] < all_tags[b]; });
    std::vector<double> sorted_xyz(3 * static_cast<std::size_t>(total));
    std::vector<int> ids(total);
    for (int k = 0; k < total; ++k) {
      const int src = order[k];
      sorted_tags[k] = all_tags[src];
      ids[k] = static_cast<int>(all_tags[src]);
      for (int d = 0; d < 3; ++d) sorted_xyz[3 * k + d] = all_xyz[3 * src + d];
    }
    const double box9[9] = {domain->boxlo[0],
                            domain->boxlo[1],
                            domain->boxlo[2],
                            domain->xprd,
                            domain->yprd,
                            domain->zprd,
                            domain->triclinic ? domain->xy : 0.0,
                            domain->triclinic ? domain->xz : 0.0,
                            domain->triclinic ? domain->yz : 0.0};
    seams_system *sys = seams_system_create();
    if (sys == nullptr) {
      status = SEAMS_ERR_NOMEM;
    } else {
      status = seams_system_set_frame(sys, sorted_xyz.data(), total, box9, nullptr, ids.data());
      if (status == SEAMS_OK) {
        status = seams_chill_plus(sys, cutoff_, 0, graph_, rule_, labels.data(), nullptr);
      }
      seams_system_destroy(sys);
    }
  }
  MPI_Bcast(&status, 1, MPI_INT, 0, world);
  if (status != SEAMS_OK) {
    error->all(FLERR, std::string("compute dseams: ") + seams_status_string(
                          static_cast<seams_status>(status)) + ": " + seams_last_error());
  }
  if (total > 0) {
    MPI_Bcast(sorted_tags.data(), total, MPI_LMP_TAGINT, 0, world);
    MPI_Bcast(labels.data(), total, MPI_INT, 0, world);
  }
  for (int i = 0; i < nlocal; ++i) {
    if (!(atom->mask[i] & groupbit)) continue;
    const auto at = std::lower_bound(sorted_tags.begin(), sorted_tags.end(), atom->tag[i]);
    label_[i] = static_cast<double>(labels[static_cast<std::size_t>(at - sorted_tags.begin())]);
  }
  seconds_ += MPI_Wtime() - t0;
  calls_ += 1;
}

double ComputeDseams::memory_usage() {
  return static_cast<double>(nmax_) * sizeof(double);
}
