/* PLUMED action that calls the d-SEAMS in-situ engine.
 *
 * Positions and the box arrive in nanometres. The engine reads angstrom.
 * SPDX-License-Identifier: MIT
 */

#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"

#include "tools/NeighborList.h"
#include "tools/Pbc.h"

#include <seams_insitu.h>

#include <memory>
#include <vector>

namespace PLMD {

class Dseams : public ActionAtomistic, public ActionPilot {
 public:
  static void registerKeywords(Keywords &keys);
  explicit Dseams(const ActionOptions &ao);
  ~Dseams() override;
  void calculate() override;
  void apply() override;

 private:
  double cutoff_;
  std::string dump_;
  std::string record_;
  seams_insitu_engine *eng_;
  std::vector<AtomNumber> atoms_;
  std::unique_ptr<NeighborList> nlist_;
};

PLUMED_REGISTER_ACTION(Dseams, "DSEAMS")

void Dseams::registerKeywords(Keywords &keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  keys.add("atoms", "ATOMS", "atoms passed to the in-situ engine");
  keys.add("compulsory", "DUMP", "dseams.lammpstrj", "trajectory written for the replay");
  keys.add("compulsory", "RECORD", "dseams.txt", "census written on each analysed step");
  keys.add("compulsory", "CUTOFF", "3.5", "neighbour cutoff in angstrom");
  keys.add("compulsory", "EVERY", "1", "analyse one frame every this many steps");
}

Dseams::Dseams(const ActionOptions &ao) :
    Action(ao), ActionAtomistic(ao), ActionPilot(ao), cutoff_(3.5), eng_(nullptr) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS", atoms);
  parse("DUMP", dump_);
  parse("RECORD", record_);
  parse("CUTOFF", cutoff_);
  int every = 1;
  parse("EVERY", every);
  checkRead();
  if (atoms.empty()) error("DSEAMS needs ATOMS");
  if (!(cutoff_ > 0.0)) error("DSEAMS cutoff must be positive");
  if (every < 1) error("DSEAMS EVERY must be positive");
  setStride(every);
  requestAtoms(atoms);
  atoms_ = atoms;
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  opt.cutoff = cutoff_;
  eng_ = seams_insitu_create(&opt);
  if (eng_ == nullptr) error("DSEAMS could not create the engine");
  if (seams_insitu_open_record(eng_, record_.c_str()) != SEAMS_INSITU_OK) error("DSEAMS could not open the record");
  log.printf("  d-SEAMS cutoff %f angstrom on %u atoms\n", cutoff_,
             static_cast<unsigned>(atoms.size()));
}

Dseams::~Dseams() { seams_insitu_destroy(eng_); }

void Dseams::apply() {}

void Dseams::calculate() {
  if (!onStep() || eng_ == nullptr) return;
  const unsigned n = getNumberOfAtoms();
  std::vector<double> xyz(3 * static_cast<std::size_t>(n));
  for (unsigned i = 0; i < n; ++i) {
    const Vector &p = getPosition(static_cast<int>(i));
    xyz[3 * i] = p[0] * 10.0;
    xyz[3 * i + 1] = p[1] * 10.0;
    xyz[3 * i + 2] = p[2] * 10.0;
  }
  const Tensor box = getPbc().getBox();
  double rows[9] = {box(0, 0) * 10.0, box(1, 0) * 10.0, box(2, 0) * 10.0,
                    box(0, 1) * 10.0, box(1, 1) * 10.0, box(2, 1) * 10.0,
                    box(0, 2) * 10.0, box(1, 2) * 10.0, box(2, 2) * 10.0};
  double xlo, ylo, zlo, lx, ly, lz, xy, xz, yz;
  if (seams_insitu_box_from_rows(rows, &xlo, &ylo, &zlo, &lx, &ly, &lz, &xy, &xz, &yz) != 0) {
    xlo = ylo = zlo = 0.0;
    lx = ly = lz = 20.0;
    xy = xz = yz = 0.0;
  }
  seams_insitu_frame frame{};
  frame.step = static_cast<int>(getStep());
  frame.n = static_cast<int>(n);
  frame.xyz = xyz.data();
  frame.xlo = xlo;
  frame.ylo = ylo;
  frame.zlo = zlo;
  frame.lx = lx;
  frame.ly = ly;
  frame.lz = lz;
  frame.xy = xy;
  frame.xz = xz;
  frame.yz = yz;
  frame.gathered = 1;
  if (nlist_ == nullptr) {
    nlist_ = std::make_unique<NeighborList>(atoms_, true, true, getPbc(), comm,
                                            (cutoff_ + 0.5) / 10.0, 1);
  }
  std::vector<Vector> pos(n);
  for (unsigned i = 0; i < n; ++i) {
    pos[i] = getPosition(static_cast<int>(i));
  }
  nlist_->update(pos);
  std::vector<int> counts(n, 1);
  std::vector<int> flat;
  flat.reserve(static_cast<std::size_t>(n) * 5);
  for (unsigned i = 0; i < n; ++i) {
    const auto neigh = nlist_->getNeighbors(i);
    counts[i] = static_cast<int>(neigh.size()) + 1;
    flat.push_back(static_cast<int>(i) + 1);
    for (const unsigned j : neigh) {
      if (j != i && j < n) {
        flat.push_back(static_cast<int>(j) + 1);
      } else {
        counts[i] -= 1;
      }
    }
  }
  seams_insitu_set_neighbours(eng_, static_cast<int>(n), counts.data(), flat.data());
  const int rc = seams_insitu_observe(eng_, &frame, nullptr);
  if (rc != SEAMS_INSITU_OK && rc != SEAMS_INSITU_SKIP) error("DSEAMS observe failed");
}

} // namespace PLMD
