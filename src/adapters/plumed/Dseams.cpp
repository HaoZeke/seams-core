/* PLUMED action that calls the d-SEAMS in-situ engine.
 *
 * Positions and the box arrive in nanometres. The engine reads angstrom.
 * SPDX-License-Identifier: MIT
 */

#include "config/version.h"
#include "core/ActionAtomistic.h"
#include "tools/Communicator.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"

#include "tools/Pbc.h"

#include <seams_insitu.h>

#include <vector>

namespace PLMD {

class Dseams : public ActionAtomistic, public ActionPilot {
 public:
  static void registerKeywords(Keywords &keys);
  explicit Dseams(const ActionOptions &ao);
  ~Dseams() override;
  void calculate() override;
  void apply() override;
#if PLUMED_VERSION_MINOR >= 10
  bool actionHasForces() override { return false; }
#endif

 private:
  double cutoff_;
  std::string dump_;
  std::string record_;
  seams_insitu_engine *eng_;
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
  if (comm.Get_rank() != 0) return;
  const unsigned n = getNumberOfAtoms();
  // Positions stay in the PLUMED array. Vector is three contiguous doubles,
  // so the array is one xyz buffer in nanometres. The engine scales it.
  const std::vector<Vector> &pos = getPositions();
  static_assert(sizeof(Vector) == 3 * sizeof(double), "PLUMED Vector is not packed xyz");
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
  frame.xyz = n > 0 ? &pos[0][0] : nullptr;
  frame.coord_scale = 10.0;
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
  // NeighborList is an all-pairs table. getIndexPair multiplies the pair
  // index by 8 in a 32-bit unsigned (NeighborList.cpp, the K formula), which
  // wraps past about 33 000 atoms and then reads off the position array.
  // The engine builds the cutoff shell with a cell list instead.
  const int rc = seams_insitu_observe(eng_, &frame, nullptr);
  if (rc != SEAMS_INSITU_OK && rc != SEAMS_INSITU_SKIP) error("DSEAMS observe failed");
}

} // namespace PLMD
