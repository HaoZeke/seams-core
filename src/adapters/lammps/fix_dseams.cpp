/* d-SEAMS LAMMPS fix.
 *
 * This file includes LAMMPS headers and is distributed under GPL-2.0-only,
 * the licence of LAMMPS. libyodaLib stays MIT and is called through
 * seams_insitu.h.
 *
 * SPDX-License-Identifier: GPL-2.0-only
 */

#include "fix.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "lammpsplugin.h"
#include "update.h"
#include "version.h"

#include <seams_insitu.h>

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace LAMMPS_NS;

namespace {

int parse_positive_int(const char *text, int *out) {
  if (text == nullptr || text[0] == '\0') return 1;
  char *end = nullptr;
  const long value = std::strtol(text, &end, 10);
  if (end == text || *end != '\0' || value < 1 || value > 1000000) {
    return 1;
  }
  *out = static_cast<int>(value);
  return 0;
}

double parse_positive_double(const char *text, int *bad) {
  if (text == nullptr || text[0] == '\0') return (*bad = 1), 0.0;
  char *end = nullptr;
  const double value = std::strtod(text, &end);
  if (end == text || *end != '\0' || !(value > 0.0)) {
    *bad = 1;
    return 0;
  }
  return value;
}

} // namespace

class FixDseams : public Fix {
 public:
  FixDseams(LAMMPS *lmp, int narg, char **arg);
  ~FixDseams() override;
  int setmask() override;
  void end_of_step() override;

 private:
  int active_;
  int every_;
  int depth_;
  double cutoff_;
  std::string dump_;
  std::string record_;
  seams_insitu_engine *eng_;
};

FixDseams::FixDseams(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), active_(0), every_(1), depth_(8), cutoff_(3.5), eng_(nullptr) {
  int iarg = 3;
  while (iarg < narg) {
    if (std::strcmp(arg[iarg], "dump") == 0 && iarg + 1 < narg) {
      dump_ = arg[iarg + 1];
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "record") == 0 && iarg + 1 < narg) {
      record_ = arg[iarg + 1];
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "every") == 0 && iarg + 1 < narg) {
      if (parse_positive_int(arg[iarg + 1], &every_) != 0) {
        active_ = 0;
        return;
      }
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "cutoff") == 0 && iarg + 1 < narg) {
      int bad = 0;
      cutoff_ = parse_positive_double(arg[iarg + 1], &bad);
      if (bad) {
        active_ = 0;
        return;
      }
      iarg += 2;
    } else if (std::strcmp(arg[iarg], "depth") == 0 && iarg + 1 < narg) {
      if (parse_positive_int(arg[iarg + 1], &depth_) != 0) {
        active_ = 0;
        return;
      }
      iarg += 2;
    } else {
      active_ = 0;
      return;
    }
  }
  if (dump_.empty() || record_.empty()) {
    active_ = 0;
    return;
  }
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.cutoff = cutoff_;
  opt.ring_depth = depth_;
  opt.cadence.every = 1;
  eng_ = seams_insitu_create(&opt);
  if (eng_ == nullptr) return;
  if (seams_insitu_open_record(eng_, record_.c_str()) != SEAMS_INSITU_OK) {
    active_ = 0;
    return;
  }
  active_ = 1;
}

FixDseams::~FixDseams() {
  seams_insitu_destroy(eng_);
  eng_ = nullptr;
}

int FixDseams::setmask() {
  int mask = 0;
  mask |= FixConst::END_OF_STEP;
  return mask;
}

void FixDseams::end_of_step() {
  if (!active_ || eng_ == nullptr) {
    return;
  }
  const auto step64 = static_cast<long long>(update->ntimestep);
  if (step64 < 0 || (step64 % every_) != 0) return;
  const int nlocal = atom->nlocal;
  std::vector<double> local(static_cast<std::size_t>(nlocal) * 3);
  std::vector<int> local_types(static_cast<std::size_t>(nlocal));
  std::vector<int> local_ids(static_cast<std::size_t>(nlocal));
  std::vector<int> mask(static_cast<std::size_t>(nlocal));
  for (int i = 0; i < nlocal; ++i) {
    local[static_cast<std::size_t>(3 * i)] = atom->x[i][0];
    local[static_cast<std::size_t>(3 * i + 1)] = atom->x[i][1];
    local[static_cast<std::size_t>(3 * i + 2)] = atom->x[i][2];
    local_types[static_cast<std::size_t>(i)] = atom->type[i];
    local_ids[static_cast<std::size_t>(i)] = static_cast<int>(atom->tag[i]);
    mask[static_cast<std::size_t>(i)] = atom->mask[i];
  }
  std::vector<double> kept_x(local.size());
  std::vector<int> kept_t(static_cast<std::size_t>(nlocal));
  std::vector<int> kept_i(static_cast<std::size_t>(nlocal));
  const int kept = seams_insitu_pack_owned(
      nlocal, mask.data(), groupbit, local.data(), local_types.data(), local_ids.data(),
      kept_x.data(), kept_t.data(), kept_i.data());
  if (kept < 0) return;
  const int nprocs = comm->nprocs;
  int counts = kept;
  int displs = 0;
  int total = 0;
  if (seams_insitu_gather_displs(1, &counts, &displs, &total) != 0 || total != kept) return;
  // One rank already holds every atom. A domain-decomposed rank sets
  // local_only and the core asks for a gather through seams_insitu_gather_displs.
  seams_insitu_frame frame{};
  frame.step = static_cast<int>(step64);
  frame.n = total;
  frame.xyz = kept_x.data();
  frame.types = kept_t.data();
  frame.ids = kept_i.data();
  frame.xlo = domain->boxlo[0];
  frame.ylo = domain->boxlo[1];
  frame.zlo = domain->boxlo[2];
  frame.lx = domain->prd[0];
  frame.ly = domain->prd[1];
  frame.lz = domain->prd[2];
  frame.xy = domain->xy;
  frame.xz = domain->xz;
  frame.yz = domain->yz;
  frame.gathered = nprocs == 1 ? 1 : 0;
  frame.local_only = nprocs > 1 ? 1 : 0;
  frame.rank = comm->me;
  frame.nranks = nprocs;
  const int rc = seams_insitu_observe(eng_, &frame, dump_.c_str());
  if (rc != SEAMS_INSITU_OK && rc != SEAMS_INSITU_SKIP) {
    active_ = 0;
  }
}

static Fix *dseams_creator(LAMMPS *lmp, int narg, char **arg) {
  return new FixDseams(lmp, narg, arg);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc) {
  auto *register_plugin = reinterpret_cast<lammpsplugin_regfunc>(regfunc);
  lammpsplugin_t plugin{};
  plugin.version = LAMMPS_VERSION;
  plugin.style = "fix";
  plugin.name = "dseams";
  plugin.info = "d-SEAMS in-situ ring census, CHILL+ and cages";
  plugin.author = "Rohit Goswami (rgoswami@ieee.org)";
  plugin.creator.v2 = reinterpret_cast<lammpsplugin_factory2 *>(&dseams_creator);
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
