// Copyright (c) d-SEAMS developers
#include "seams_lang.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <sys/resource.h>
#include <unistd.h>

namespace {

int fail(const std::string &message) {
  std::cerr << message;
  const char *err = seams_last_error();
  if (err != nullptr && err[0] != '\0') {
    std::cerr << ": " << err;
  }
  std::cerr << "\n";
  return 1;
}

// Every driver check calls this. A false condition is the failure arm.
// The self-test below forces that arm and checks the text it prints.
int check(bool ok, const char *message) {
  if (!ok) {
    return fail(message);
  }
  return 0;
}


int probe_arm(const std::function<int(bool)> &arm, const char *message) {
  int pipefd[2];
  if (pipe(pipefd) != 0) {
    return 1;
  }
  const int saved = dup(STDERR_FILENO);
  dup2(pipefd[1], STDERR_FILENO);
  const int status = arm(true);
  std::cerr.flush();
  dup2(saved, STDERR_FILENO);
  close(saved);
  close(pipefd[1]);
  std::string text;
  char buf[512];
  ssize_t nread = 0;
  while ((nread = read(pipefd[0], buf, sizeof buf)) > 0) {
    text.append(buf, static_cast<std::size_t>(nread));
  }
  close(pipefd[0]);
  if (status != 1) {
    return 1;
  }
  if (text.find(message) == std::string::npos) {
    return 1;
  }
  return 0;
}

int self_test_fail_arm() {
  int pipefd[2];
  if (pipe(pipefd) != 0) {
    return check(false, "fail arm pipe");
  }
  const int saved = dup(STDERR_FILENO);
  dup2(pipefd[1], STDERR_FILENO);
  (void)seams_path_nframes(nullptr);
  const int status = check(false, "fail arm");
  std::cerr.flush();
  dup2(saved, STDERR_FILENO);
  close(saved);
  close(pipefd[1]);
  std::string text;
  char buf[256];
  ssize_t nread = 0;
  while ((nread = read(pipefd[0], buf, sizeof buf)) > 0) {
    text.append(buf, static_cast<std::size_t>(nread));
  }
  close(pipefd[0]);
  int bad = 0;
    {
  auto arm0 = [&](bool force) -> int {
    if (force || (!(status == 1))) {
      return fail("fail arm status");
    }
    return 0;
  };
  bad += probe_arm(arm0, "fail arm status");
  bad += arm0(false) != 0;
  }
    {
  auto arm1 = [&](bool force) -> int {
    if (force || (!(text.find("fail arm") != std::string::npos))) {
      return fail("fail arm text");
    }
    return 0;
  };
  bad += probe_arm(arm1, "fail arm text");
  bad += arm1(false) != 0;
  }
    {
  auto arm2 = [&](bool force) -> int {
    if (force || (!(text.find("path is null") != std::string::npos))) {
      return fail("fail arm error");
    }
    return 0;
  };
  bad += probe_arm(arm2, "fail arm error");
  bad += arm2(false) != 0;
  }
  if (bad != 0) {
    return 1;
  }
  return 0;
}

int count_of(const int counts[9], int label) { return counts[label]; }

} // namespace

int main(int argc, char **argv) {
  const std::string root = argc > 1 ? argv[1] : ".";
  const std::string tip4p = root + "/input/traj/exampleTraj.lammpstrj";
  const std::string cubic = root + "/input/traj/mW_cubic.lammpstrj";
  const std::string si = root + "/input/traj/genice_sI.lammpstrj";
  if (self_test_fail_arm() != 0) {
    return 1;
  }
  int bad = 0;
  {
 auto arm3 = [&](bool force) -> int {
   if (force || (std::string(seams_lang_version()).empty())) {
     return fail("empty version");
   }
   return 0;
 };
 bad += probe_arm(arm3, "empty version");
 bad += arm3(false) != 0;
 }

  seams_frame *tip = seams_frame_read(tip4p.c_str(), 1, 2);
    {
  auto arm4 = [&](bool force) -> int {
    if (force || (tip == nullptr)) {
      return fail("read tip4p");
    }
    return 0;
  };
  bad += probe_arm(arm4, "read tip4p");
  bad += arm4(false) != 0;
  }
    {
  auto arm5 = [&](bool force) -> int {
    if (force || (seams_frame_nop(tip) != 250)) {
      return fail("tip4p nop");
    }
    return 0;
  };
  bad += probe_arm(arm5, "tip4p nop");
  bad += arm5(false) != 0;
  }
  const double *xyz = seams_frame_positions(tip);
    {
  auto arm6 = [&](bool force) -> int {
    if (force || (xyz == nullptr)) {
      return fail("positions");
    }
    return 0;
  };
  bad += probe_arm(arm6, "positions");
  bad += arm6(false) != 0;
  }

  std::vector<int> labels(250);
  int counts[9] = {};
    {
  auto arm7 = [&](bool force) -> int {
    if (force || (seams_frame_chill_plus(tip, 3.5, 2, "mutual", labels.data(), counts) != 0)) {
      return fail("chill+ mutual");
    }
    return 0;
  };
  bad += probe_arm(arm7, "chill+ mutual");
  bad += arm7(false) != 0;
  }
    {
  auto arm8 = [&](bool force) -> int {
    if (force || (count_of(counts, 4) != 1 || count_of(counts, 5) != 7 || count_of(counts, 2) != 242)) {
      return fail("mutual chill+ counts");
    }
    return 0;
  };
  bad += probe_arm(arm8, "mutual chill+ counts");
  bad += arm8(false) != 0;
  }

    {
  auto arm9 = [&](bool force) -> int {
    if (force || (seams_frame_chill_plus(tip, 3.5, 2, "cutoff", labels.data(), counts) != 0)) {
      return fail("chill+ cutoff");
    }
    return 0;
  };
  bad += probe_arm(arm9, "chill+ cutoff");
  bad += arm9(false) != 0;
  }
    {
  auto arm10 = [&](bool force) -> int {
    if (force || (count_of(counts, 5) != 12 || count_of(counts, 2) != 238 || count_of(counts, 4) != 0)) {
      return fail("cutoff chill+ counts");
    }
    return 0;
  };
  bad += probe_arm(arm10, "cutoff chill+ counts");
  bad += arm10(false) != 0;
  }
  int default_counts[9] = {};
    {
  auto arm11 = [&](bool force) -> int {
    if (force || (seams_frame_chill_plus(tip, 3.5, 2, nullptr, labels.data(), default_counts) != 0)) {
      return fail("chill+ default");
    }
    return 0;
  };
  bad += probe_arm(arm11, "chill+ default");
  bad += arm11(false) != 0;
  }
  for (int i = 0; i < 9; ++i) {
        {
    auto arm12 = [&](bool force) -> int {
      if (force || (default_counts[i] != counts[i])) {
        return fail("chill+ default is not the cutoff list");
      }
      return 0;
    };
    bad += probe_arm(arm12, "chill+ default is not the cutoff list");
    bad += arm12(false) != 0;
    }
  }
    {
  auto arm13 = [&](bool force) -> int {
    if (force || (seams_frame_chill_plus(tip, 3.5, 2, "knn-union", labels.data(), counts) == 0)) {
      return fail("chill+ accepted knn-union");
    }
    return 0;
  };
  bad += probe_arm(arm13, "chill+ accepted knn-union");
  bad += arm13(false) != 0;
  }

  int n_hex = -1, n_cubic = -1, n_water = -1, n_both = -1;
    {
  auto arm14 = [&](bool force) -> int {
    if (force || (seams_frame_cages(tip, 3.5, 2, 4, "seeded", 0, nullptr, &n_hex, &n_cubic, &n_water, &n_both) != 0)) {
      return fail("cages");
    }
    return 0;
  };
  bad += probe_arm(arm14, "cages");
  bad += arm14(false) != 0;
  }
    {
  auto arm15 = [&](bool force) -> int {
    if (force || (n_hex != 0 || n_cubic != 0 || n_water != 250)) {
      return fail("cage counts");
    }
    return 0;
  };
  bad += probe_arm(arm15, "cage counts");
  bad += arm15(false) != 0;
  }

  std::vector<double> ql(250), qlbar(250);
    {
  auto arm16 = [&](bool force) -> int {
    if (force || (seams_frame_steinhardt(tip, 3.5, 2, 3, ql.data(), qlbar.data()) != 0)) {
      return fail("steinhardt");
    }
    return 0;
  };
  bad += probe_arm(arm16, "steinhardt");
  bad += arm16(false) != 0;
  }
  double mean = 0.0;
  int nfinite = 0;
  for (double value : ql) {
    if (std::isfinite(value)) {
      mean += value;
      ++nfinite;
    }
  }
  mean = nfinite > 0 ? mean / static_cast<double>(nfinite) : 0.0;
    {
  auto arm17 = [&](bool force) -> int {
    if (force || (std::abs(mean - 0.389362) > 1e-4)) {
      return fail("q3");
    }
    return 0;
  };
  bad += probe_arm(arm17, "q3");
  bad += arm17(false) != 0;
  }

  const int bins = 4;
  std::vector<double> r(bins), g(bins);
  std::vector<int> hist(bins);
  double volume = 0.0;
    {
  auto arm18 = [&](bool force) -> int {
    if (force || (seams_frame_rdf(tip, 6.0, bins, 2, 2, r.data(), g.data(), hist.data(), &volume) != 0)) {
      return fail("rdf");
    }
    return 0;
  };
  bad += probe_arm(arm18, "rdf");
  bad += arm18(false) != 0;
  }
    {
  auto arm19 = [&](bool force) -> int {
    if (force || (hist[1] != 390 || hist[2] != 705 || hist[3] != 531)) {
      return fail("rdf counts");
    }
    return 0;
  };
  bad += probe_arm(arm19, "rdf counts");
  bad += arm19(false) != 0;
  }

  std::vector<int> rings(9, 0);
    {
  auto arm20 = [&](bool force) -> int {
    if (force || (seams_frame_rings(tip, 3.5, 2, "mutual", 8, rings.data(), 9) != 0)) {
      return fail("rings");
    }
    return 0;
  };
  bad += probe_arm(arm20, "rings");
  bad += arm20(false) != 0;
  }
  std::cout << "tip4p mutual rings";
  for (int i = 3; i <= 8; ++i) {
    std::cout << " " << i << ":" << rings[i];
  }
  std::cout << "\n";

  seams_frame *ice = seams_frame_read(cubic.c_str(), 1, 1);
    {
  auto arm21 = [&](bool force) -> int {
    if (force || (ice == nullptr || seams_frame_nop(ice) != 4096)) {
      return fail("mW");
    }
    return 0;
  };
  bad += probe_arm(arm21, "mW");
  bad += arm21(false) != 0;
  }
  std::vector<int> ice_labels(4096);
  int ice_counts[9] = {};
    {
  auto arm22 = [&](bool force) -> int {
    if (force || (seams_frame_chill_plus(ice, 3.5, 1, "mutual", ice_labels.data(), ice_counts) != 0)) {
      return fail("mW chill+");
    }
    return 0;
  };
  bad += probe_arm(arm22, "mW chill+");
  bad += arm22(false) != 0;
  }
    {
  auto arm23 = [&](bool force) -> int {
    if (force || (ice_counts[0] != 4096)) {
      return fail("mW cubic count");
    }
    return 0;
  };
  bad += probe_arm(arm23, "mW cubic count");
  bad += arm23(false) != 0;
  }
    {
  auto arm24 = [&](bool force) -> int {
    if (force || (seams_frame_cages(ice, 3.5, 1, 4, "seeded", 0, nullptr, &n_hex, &n_cubic, &n_water, &n_both) != 0)) {
      return fail("mW cages");
    }
    return 0;
  };
  bad += probe_arm(arm24, "mW cages");
  bad += arm24(false) != 0;
  }
    {
  auto arm25 = [&](bool force) -> int {
    if (force || (n_hex != 0 || n_cubic != 4096 || n_water != 0)) {
      return fail("mW cage counts");
    }
    return 0;
  };
  bad += probe_arm(arm25, "mW cage counts");
  bad += arm25(false) != 0;
  }
    {
  auto arm26 = [&](bool force) -> int {
    if (force || (seams_path_nframes(cubic.c_str()) < 2)) {
      return fail("mW nframes");
    }
    return 0;
  };
  bad += probe_arm(arm26, "mW nframes");
  bad += arm26(false) != 0;
  }

  seams_frame *guest = seams_frame_read(si.c_str(), 1, 1);
    {
  auto arm27 = [&](bool force) -> int {
    if (force || (guest == nullptr)) {
      return fail("sI");
    }
    return 0;
  };
  bad += probe_arm(arm27, "sI");
  bad += arm27(false) != 0;
  }
  int n_cages = 0;
  int n_atoms = 0;
    {
  auto arm28 = [&](bool force) -> int {
    if (force || (seams_frame_signature(guest, 3.5, 1, "cutoff", "512", &n_cages, &n_atoms) != 0)) {
      return fail("sI signature");
    }
    return 0;
  };
  bad += probe_arm(arm28, "sI signature");
  bad += arm28(false) != 0;
  }
    {
  auto arm29 = [&](bool force) -> int {
    if (force || (seams_frame_nop(guest) != 46 || n_cages != 2 || n_atoms != 40)) {
      return fail("sI 512 census");
    }
    return 0;
  };
  bad += probe_arm(arm29, "sI 512 census");
  bad += arm29(false) != 0;
  }
    {
  auto arm30 = [&](bool force) -> int {
    if (force || (seams_frame_signature_edge(guest, 3.5, 1, "cutoff", "512", "none",
                                 &n_cages, &n_atoms) != 0 || n_cages != 2 || n_atoms != 40)) {
      return fail("sI 512 none");
    }
    return 0;
  };
  bad += probe_arm(arm30, "sI 512 none");
  bad += arm30(false) != 0;
  }
    {
  auto arm31 = [&](bool force) -> int {
    if (force || (seams_frame_signature_edge(guest, 3.5, 1, "cutoff", "512", "1,2",
                                 &n_cages, &n_atoms) != 0 || n_cages != 0)) {
      return fail("sI 512 tight window");
    }
    return 0;
  };
  bad += probe_arm(arm31, "sI 512 tight window");
  bad += arm31(false) != 0;
  }
    {
  auto arm32 = [&](bool force) -> int {
    if (force || (seams_frame_signature(guest, 3.5, 1, "cutoff", "51268", &n_cages, &n_atoms) != 0)) {
      return fail("sI 51268");
    }
    return 0;
  };
  bad += probe_arm(arm32, "sI 51268");
  bad += arm32(false) != 0;
  }
    {
  auto arm33 = [&](bool force) -> int {
    if (force || (n_cages != 0)) {
      return fail("sI 51268 census");
    }
    return 0;
  };
  bad += probe_arm(arm33, "sI 51268 census");
  bad += arm33(false) != 0;
  }

  const std::string nucleation = root + "/input/traj/nucleation_f1.lammpstrj";
  seams_frame *molten = seams_frame_read(nucleation.c_str(), 1, 1);
    {
  auto arm34 = [&](bool force) -> int {
    if (force || (molten == nullptr)) {
      return fail("nucleation frame");
    }
    return 0;
  };
  bad += probe_arm(arm34, "nucleation frame");
  bad += arm34(false) != 0;
  }
  const int nMolten = seams_frame_nop(molten);
  std::vector<int> hostLabels(static_cast<std::size_t>(nMolten));
  std::vector<int> offLabels(static_cast<std::size_t>(nMolten));
  int hHex = 0, hCub = 0, hWat = 0, hBoth = 0;
  int oHex = 0, oCub = 0, oWat = 0, oBoth = 0;
  unsetenv("SEAMS_OFFLOAD");
    {
  auto arm35 = [&](bool force) -> int {
    if (force || (seams_frame_cages(molten, 3.5, 1, 4, "seeded", 0, hostLabels.data(), &hHex, &hCub, &hWat, &hBoth) != 0)) {
      return fail("seeded host");
    }
    return 0;
  };
  bad += probe_arm(arm35, "seeded host");
  bad += arm35(false) != 0;
  }
  setenv("SEAMS_OFFLOAD", "1", 1);
    {
  auto arm36 = [&](bool force) -> int {
    if (force || (seams_frame_cages(molten, 3.5, 1, 4, "seeded", 0, offLabels.data(), &oHex, &oCub, &oWat, &oBoth) != 0)) {
      return fail("seeded offload");
    }
    return 0;
  };
  bad += probe_arm(arm36, "seeded offload");
  bad += arm36(false) != 0;
  }
  unsetenv("SEAMS_OFFLOAD");
    {
  auto arm37 = [&](bool force) -> int {
    if (force || (hHex != oHex || hCub != oCub || hWat != oWat || hBoth != oBoth || hostLabels != offLabels)) {
      return fail("seeded offload disagrees with the host rule");
    }
    return 0;
  };
  bad += probe_arm(arm37, "seeded offload disagrees with the host rule");
  bad += arm37(false) != 0;
  }
    {
  auto arm38 = [&](bool force) -> int {
    if (force || (hHex + hCub + hBoth + hWat != nMolten)) {
      return fail("disjoint cage counters");
    }
    return 0;
  };
  bad += probe_arm(arm38, "disjoint cage counters");
  bad += arm38(false) != 0;
  }
  int codes[4] = {};
  for (int code : hostLabels) {
        {
    auto arm39 = [&](bool force) -> int {
      if (force || (code < 0 || code > 3)) {
        return fail("cage code");
      }
      return 0;
    };
    bad += probe_arm(arm39, "cage code");
    bad += arm39(false) != 0;
    }
    ++codes[code];
  }
  // 0 water, 1 hexagonal only, 2 cubic only, 3 both
    {
  auto arm40 = [&](bool force) -> int {
    if (force || (codes[1] != hHex || codes[2] != hCub || codes[3] != hBoth || codes[0] != hWat || hBoth <= 0)) {
      return fail("stacking-fault counters");
    }
    return 0;
  };
  bad += probe_arm(arm40, "stacking-fault counters");
  bad += arm40(false) != 0;
  }
  const std::string cli = argc > 2 ? std::string(argv[2]) : std::string("bbdir-cov/seams");
  if (!cli.empty()) {
    const std::string cmd = cli + " cages " + nucleation +
                            " --graph seeded --cutoff 3.5";
    FILE *pipe = popen(cmd.c_str(), "r");
        {
    auto arm41 = [&](bool force) -> int {
      if (force || (pipe == nullptr)) {
        return fail("cli pipe");
      }
      return 0;
    };
    bad += probe_arm(arm41, "cli pipe");
    bad += arm41(false) != 0;
    }
    std::string line;
    char buf[512];
    while (fgets(buf, sizeof buf, pipe) != nullptr) {
      line += buf;
    }
    const int rc = pclose(pipe);
        {
    auto arm42 = [&](bool force) -> int {
      if (force || (rc != 0)) {
        return fail("cli cages");
      }
      return 0;
    };
    bad += probe_arm(arm42, "cli cages");
    bad += arm42(false) != 0;
    }
    int cliNop = 0, cliHex = 0, cliCub = 0, cliBoth = 0, cliWat = 0;
    std::istringstream in(line);
    std::string tok;
    in >> tok >> cliNop >> tok >> tok >> tok >> cliHex >> tok >> cliCub >>
        tok >> cliBoth >> tok >> cliWat;
        {
    auto arm43 = [&](bool force) -> int {
      if (force || (cliHex != hHex || cliCub != hCub || cliBoth != hBoth || cliWat != hWat)) {
        return fail("abi disagrees with the cli");
      }
      return 0;
    };
    bad += probe_arm(arm43, "abi disagrees with the cli");
    bad += arm43(false) != 0;
    }
  }
  {
 auto arm44 = [&](bool force) -> int {
   if (force || (seams_last_error() == nullptr)) {
     return fail("error pointer");
   }
   return 0;
 };
 bad += probe_arm(arm44, "error pointer");
 bad += arm44(false) != 0;
 }
    {
  auto arm45 = [&](bool force) -> int {
    if (force || (seams_path_nframes(nullptr) != -1 || seams_frame_read(nullptr, 1, -1) != nullptr)) {
      return fail("null path");
    }
    return 0;
  };
  bad += probe_arm(arm45, "null path");
  bad += arm45(false) != 0;
  }
  double box9[9] = {};
  seams_frame_box(nullptr, box9);
  seams_frame_box(tip, nullptr);
  seams_frame_box(tip, box9);
    {
  auto arm46 = [&](bool force) -> int {
    if (force || (seams_frame_index(nullptr) != 0 || seams_frame_nop(nullptr) != 0)) {
      return fail("null frame queries");
    }
    return 0;
  };
  bad += probe_arm(arm46, "null frame queries");
  bad += arm46(false) != 0;
  }
    {
  auto arm47 = [&](bool force) -> int {
    if (force || (seams_frame_index(tip) < 1 || seams_frame_types(tip) == nullptr ||
      seams_frame_ids(tip) == nullptr || seams_frame_types(nullptr) != nullptr || seams_frame_ids(nullptr) != nullptr || seams_frame_positions(nullptr) != nullptr)) {
      return fail("frame columns");
    }
    return 0;
  };
  bad += probe_arm(arm47, "frame columns");
  bad += arm47(false) != 0;
  }
    {
  auto arm48 = [&](bool force) -> int {
    if (force || (seams_frame_from_arrays(nullptr, 4, box9, nullptr, nullptr) != nullptr)) {
      return fail("null arrays");
    }
    return 0;
  };
  bad += probe_arm(arm48, "null arrays");
  bad += arm48(false) != 0;
  }
  const double builtXyz[12] = {0, 0, 0, 1.2, 0, 0, 0, 1.2, 0, 0.4, 0.4, 1.2};
  const double ortho[9] = {0, 0, 0, 8, 8, 8, 0, 0, 0};
  const double sheared[9] = {0, 0, 0, 8, 8, 8, 1.5, 0, 0};
  const int types4[4] = {1, 1, 1, 1};
  seams_frame *built = seams_frame_from_arrays(builtXyz, 4, ortho, types4, nullptr);
  seams_frame *tilted = seams_frame_from_arrays(builtXyz, 4, sheared, nullptr, nullptr);
    {
  auto arm49 = [&](bool force) -> int {
    if (force || (built == nullptr || tilted == nullptr)) {
      return fail("from arrays");
    }
    return 0;
  };
  bad += probe_arm(arm49, "from arrays");
  bad += arm49(false) != 0;
  }
  int labels4[4] = {};
  int counts9[9] = {};
    {
  auto arm50 = [&](bool force) -> int {
    if (force || (seams_frame_chill(built, 3.5, 1, "cutoff", labels4, counts9) != 0 ||
      seams_frame_chill(nullptr, 3.5, 1, "mutual", labels4, counts9) == 0 || seams_frame_chill(built, 3.5, 0, "nope", labels4, counts9) == 0)) {
      return fail("chill");
    }
    return 0;
  };
  bad += probe_arm(arm50, "chill");
  bad += arm50(false) != 0;
  }
  int ringCounts[9] = {};
    {
  auto arm51 = [&](bool force) -> int {
    if (force || (seams_frame_rings(built, 3.0, 1, "knn-union", 6, ringCounts, 9) != 0 ||
      seams_frame_rings(built, 3.0, 1, "cutoff", 6, ringCounts, 9) != 0 ||
      seams_frame_rings(built, 3.0, 1, "knn", 6, ringCounts, 9) != 0 ||
      seams_frame_rings(built, 3.0, 0, nullptr, 6, ringCounts, 9) != 0 ||
      seams_frame_rings(nullptr, 3.0, 1, "mutual", 6, ringCounts, 9) == 0 ||
      seams_frame_rings(built, 3.0, 1, "seeded", 6, ringCounts, 9) == 0 || seams_frame_rings(built, 3.0, 1, "nope", 6, ringCounts, 9) == 0)) {
      return fail("rings edges");
    }
    return 0;
  };
  bad += probe_arm(arm51, "rings edges");
  bad += arm51(false) != 0;
  }
    {
  auto arm52 = [&](bool force) -> int {
    if (force || (seams_frame_rdf_nbins(6.0, 0) < 1 || seams_frame_rdf2d_nbin(3.0, 0.0) < 1 || seams_frame_rdf2d_nbin(0.01, 1.0) < 1)) {
      return fail("bin counts");
    }
    return 0;
  };
  bad += probe_arm(arm52, "bin counts");
  bad += arm52(false) != 0;
  }
  // nbins 0 asks for rmax/0.1 bins. 4.0 / 0.1 is 40.
  double radius[64] = {};
  double gram[64] = {};
  int histEdge[64] = {};
  double volumeEdge = 0.0;
    {
  auto arm53 = [&](bool force) -> int {
    if (force || (seams_frame_rdf(built, 4.0, 0, 1, 0, radius, gram, histEdge, &volumeEdge) != 0 || seams_frame_rdf(nullptr, 4.0, 4, 1, 1, radius, gram, histEdge, nullptr) == 0)) {
      return fail("rdf edges");
    }
    return 0;
  };
  bad += probe_arm(arm53, "rdf edges");
  bad += arm53(false) != 0;
  }
  const std::string outdir = "var/run/dseams-lang-out";
    {
  auto arm54 = [&](bool force) -> int {
    if (force || (seams_frame_rdf2d(built, 3.0, 0.5, nullptr, radius, gram) != 0 ||
      seams_frame_rdf2d(built, 3.0, 0.5, outdir.c_str(), radius, gram) != 0 || seams_frame_rdf2d(nullptr, 3.0, 0.5, nullptr, radius, gram) == 0)) {
      return fail("rdf2d");
    }
    return 0;
  };
  bad += probe_arm(arm54, "rdf2d");
  bad += arm54(false) != 0;
  }
  int nrows = 0;
    {
  auto arm55 = [&](bool force) -> int {
    if (force || (seams_frame_prisms(built, 3.0, 1, 6, nullptr, nullptr, &nrows) != 0 || nrows < 1)) {
      return fail("prism count");
    }
    return 0;
  };
  bad += probe_arm(arm55, "prism count");
  bad += arm55(false) != 0;
  }
  std::vector<seams_prism_row> prows(static_cast<std::size_t>(nrows));
    {
  auto arm56 = [&](bool force) -> int {
    if (force || (seams_frame_prisms(nullptr, 3.0, 1, 2, nullptr, nullptr, nullptr) == 0 || seams_frame_prisms(built, 3.0, 1, 6, outdir.c_str(), prows.data(), &nrows) != 0)) {
      return fail("prisms");
    }
    return 0;
  };
  bad += probe_arm(arm56, "prisms");
  bad += arm56(false) != 0;
  }
  int mono[9] = {};
    {
  auto arm57 = [&](bool force) -> int {
    if (force || (seams_frame_monolayer(built, 3.0, 1, 6, 10.0, nullptr, mono, 9) != 0 ||
      seams_frame_monolayer(built, 3.0, 1, 6, 10.0, outdir.c_str(), mono, 9) != 0 || seams_frame_monolayer(nullptr, 3.0, 1, 6, 10.0, nullptr, mono, 9) == 0)) {
      return fail("monolayer");
    }
    return 0;
  };
  bad += probe_arm(arm57, "monolayer");
  bad += arm57(false) != 0;
  }
  int nHex = 0, nCub = 0, nWat = 0, nBoth = 0;
  int perAtom[4] = {};
    {
  auto arm58 = [&](bool force) -> int {
    if (force || (seams_frame_cages(built, 3.0, 1, 4, "cutoff", 1, perAtom, &nHex, &nCub, &nWat, &nBoth) != 0 || seams_frame_cages(nullptr, 3.0, 1, 4, "seeded", 0, nullptr, &nHex, &nCub, &nWat, &nBoth) == 0)) {
      return fail("cage edges");
    }
    return 0;
  };
  bad += probe_arm(arm58, "cage edges");
  bad += arm58(false) != 0;
  }
  double qlEdge[4] = {};
  double qlbarEdge[4] = {};
    {
  auto arm59 = [&](bool force) -> int {
    if (force || (seams_frame_steinhardt(nullptr, 3.0, 1, 4, qlEdge, qlbarEdge) == 0 || seams_frame_steinhardt(built, 3.0, 1, 4, nullptr, nullptr) != 0)) {
      return fail("steinhardt edges");
    }
    return 0;
  };
  bad += probe_arm(arm59, "steinhardt edges");
  bad += arm59(false) != 0;
  }
  int nCages = 0, nAtoms = 0;
    {
  auto arm60 = [&](bool force) -> int {
    if (force || (seams_frame_signature(nullptr, 3.0, 1, "cutoff", "512", &nCages, &nAtoms) == 0 || seams_frame_signature(built, 3.0, 1, "nope", "512", &nCages, &nAtoms) == 0)) {
      return fail("signature edges");
    }
    return 0;
  };
  bad += probe_arm(arm60, "signature edges");
  bad += arm60(false) != 0;
  }
  double tiltedBox[9] = {};
  seams_frame_box(tilted, tiltedBox);
    {
  auto arm61 = [&](bool force) -> int {
    if (force || (tiltedBox[6] == 0.0)) {
      return fail("shear was dropped");
    }
    return 0;
  };
  bad += probe_arm(arm61, "shear was dropped");
  bad += arm61(false) != 0;
  }
    {
  auto arm62 = [&](bool force) -> int {
    if (force || (seams_path_nframes("nodot") < 0 && seams_last_error() == nullptr)) {
      return fail("extension");
    }
    return 0;
  };
  bad += probe_arm(arm62, "extension");
  bad += arm62(false) != 0;
  }
  int guestHex = 0, guestCub = 0, guestWat = 0, guestBoth = 0;
    {
  auto arm63 = [&](bool force) -> int {
    if (force || (seams_frame_cages(guest, 3.5, 1, 0, "cutoff", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) != 0)) {
      return fail("guest cutoff cages");
    }
    return 0;
  };
  bad += probe_arm(arm63, "guest cutoff cages");
  bad += arm63(false) != 0;
  }
  setenv("SEAMS_OFFLOAD", "1", 1);
    {
  auto arm64 = [&](bool force) -> int {
    if (force || (seams_frame_cages(guest, 3.5, 1, 4, "mutual", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) != 0)) {
      return fail("guest offload cages");
    }
    return 0;
  };
  bad += probe_arm(arm64, "guest offload cages");
  bad += arm64(false) != 0;
  }
  unsetenv("SEAMS_OFFLOAD");
  const std::string emptyDump = (std::filesystem::temp_directory_path() / "dseams-lang-empty.lammpstrj").string();
  {
    std::ofstream out(emptyDump);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
           "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
           "ITEM: ATOMS id type x y z\n";
  }
  seams_frame *empty = seams_frame_read(emptyDump.c_str(), 1, -1);
    {
  auto arm65 = [&](bool force) -> int {
    if (force || (empty == nullptr || seams_frame_cages(empty, 3.0, 0, 4, "cutoff", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) != 0)) {
      return fail("empty cages");
    }
    return 0;
  };
  bad += probe_arm(arm65, "empty cages");
  bad += arm65(false) != 0;
  }
  const std::string tinyXyz = (std::filesystem::temp_directory_path() / "dseams-lang-tiny.xyz").string();
  {
    std::ofstream out(tinyXyz);
    out << "1\ncomment\nO 0 0 0\n";
  }
  seams_frame *fromXyz = seams_frame_read(tinyXyz.c_str(), 1, -1);
    {
  auto arm66 = [&](bool force) -> int {
    if (force || (fromXyz == nullptr)) {
      return fail("xyz read");
    }
    return 0;
  };
  bad += probe_arm(arm66, "xyz read");
  bad += arm66(false) != 0;
  }
  const auto nroot = std::filesystem::temp_directory_path();
  const std::string twoXyz = (nroot / "dseams-nframes-two.xyz").string();
  const std::string oneGro = (nroot / "dseams-nframes.gro").string();
  const std::string onePdb = (nroot / "dseams-nframes.pdb").string();
  {
    std::ofstream out(twoXyz);
    out << "1\nframe\nO 0.0 0.0 0.0\n1\nframe\nO 1.0 0.0 0.0\n";
  }
  {
    std::ofstream out(oneGro);
    out << "water\n    1\n"
           "    1SOL     OW    1   0.000   0.000   0.000\n"
           "   1.00000   1.00000   1.00000\n";
  }
  {
    std::ofstream out(onePdb);
    out << "ATOM      1  OW  SOL     1       0.000   0.000   0.000  1.00  0.00           O\nEND\n";
  }
    {
  auto arm67 = [&](bool force) -> int {
    if (force || (seams_path_nframes(twoXyz.c_str()) != 2)) {
      return fail("xyz frame count");
    }
    return 0;
  };
  bad += probe_arm(arm67, "xyz frame count");
  bad += arm67(false) != 0;
  }
    {
  auto arm68 = [&](bool force) -> int {
    if (force || (seams_path_nframes(oneGro.c_str()) != 1)) {
      return fail("gro frame count");
    }
    return 0;
  };
  bad += probe_arm(arm68, "gro frame count");
  bad += arm68(false) != 0;
  }
    {
  auto arm69 = [&](bool force) -> int {
    if (force || (seams_path_nframes(onePdb.c_str()) != 1)) {
      return fail("pdb frame count");
    }
    return 0;
  };
  bad += probe_arm(arm69, "pdb frame count");
  bad += arm69(false) != 0;
  }
  const std::string missingXyzCount = (nroot / "dseams-nframes-missing.xyz").string();
  const std::string badXyz = (nroot / "dseams-nframes-bad.xyz").string();
  {
    std::ofstream out(badXyz);
    out << "this is not an xyz file\n";
  }
    {
  auto arm70 = [&](bool force) -> int {
    if (force || (seams_path_nframes(missingXyzCount.c_str()) != -1 || seams_last_error() == nullptr)) {
      return fail("missing xyz count");
    }
    return 0;
  };
  bad += probe_arm(arm70, "missing xyz count");
  bad += arm70(false) != 0;
  }
    {
  auto arm71 = [&](bool force) -> int {
    if (force || (seams_path_nframes(badXyz.c_str()) != -1 || seams_last_error() == nullptr)) {
      return fail("bad xyz count");
    }
    return 0;
  };
  bad += probe_arm(arm71, "bad xyz count");
  bad += arm71(false) != 0;
  }
    {
  auto arm72 = [&](bool force) -> int {
    if (force || (seams_path_nframes("missing-file.lammpstrj") != -1 || seams_last_error() == nullptr)) {
      return fail("missing dump count");
    }
    return 0;
  };
  bad += probe_arm(arm72, "missing dump count");
  bad += arm72(false) != 0;
  }
  // A thrown error must come back as a status, not abort the host.
  auto tighten = []() {
    std::ifstream status("/proc/self/status");
    std::string key;
    long kb = 0;
    while (status >> key) {
      if (key == "VmSize:") {
        status >> kb;
        break;
      }
      status.ignore(1 << 20, '\n');
    }
    struct rlimit lim {};
    getrlimit(RLIMIT_AS, &lim);
    const rlim_t old = lim.rlim_cur;
    lim.rlim_cur = static_cast<rlim_t>(kb) * 1024 + 8 * 1024 * 1024;
    if (setrlimit(RLIMIT_AS, &lim) != 0) {
      return old;
    }
    return old;
  };
  auto loosen = [](rlim_t old) {
    struct rlimit lim {};
    getrlimit(RLIMIT_AS, &lim);
    lim.rlim_cur = old;
    setrlimit(RLIMIT_AS, &lim);
  };
  const int many = 2000000;
  std::vector<double> wide(static_cast<std::size_t>(many) * 3, 0.0);
  double wideBox[9] = {0, 0, 0, 10, 10, 10, 0, 0, 0};
  rlim_t cap = tighten();
  seams_frame *over = seams_frame_from_arrays(wide.data(), many, wideBox, nullptr, nullptr);
  loosen(cap);
    {
  auto arm73 = [&](bool force) -> int {
    if (force || (over != nullptr || seams_last_error() == nullptr)) {
      return fail("array guard");
    }
    return 0;
  };
  bad += probe_arm(arm73, "array guard");
  bad += arm73(false) != 0;
  }
  seams_frame_free(over);
  std::vector<double> rr(static_cast<std::size_t>(many));
  std::vector<double> gg(static_cast<std::size_t>(many));
  std::vector<int> cc(static_cast<std::size_t>(many));
  cap = tighten();
  const int rdfStatus = seams_frame_rdf(built, 10.0, many, 1, 1, rr.data(), gg.data(), cc.data(), nullptr);
  loosen(cap);
    {
  auto arm74 = [&](bool force) -> int {
    if (force || (rdfStatus == 0 || seams_last_error() == nullptr)) {
      return fail("rdf guard");
    }
    return 0;
  };
  bad += probe_arm(arm74, "rdf guard");
  bad += arm74(false) != 0;
  }
  cap = tighten();
  const int rdf2Status = seams_frame_rdf2d(built, 1.0e8, 0.05, nullptr, rr.data(), gg.data());
  loosen(cap);
    {
  auto arm75 = [&](bool force) -> int {
    if (force || (rdf2Status == 0 || seams_last_error() == nullptr)) {
      return fail("rdf2d guard");
    }
    return 0;
  };
  bad += probe_arm(arm75, "rdf2d guard");
  bad += arm75(false) != 0;
  }
  seams_frame *crowd = seams_frame_from_arrays(wide.data(), many, wideBox, nullptr, nullptr);
    {
  auto arm76 = [&](bool force) -> int {
    if (force || (crowd == nullptr)) {
      return fail("crowd frame");
    }
    return 0;
  };
  bad += probe_arm(arm76, "crowd frame");
  bad += arm76(false) != 0;
  }
  cap = tighten();
  const int steinStatus = seams_frame_steinhardt(crowd, 3.5, 1, 4, rr.data(), gg.data());
  loosen(cap);
  seams_frame_free(crowd);
    {
  auto arm77 = [&](bool force) -> int {
    if (force || (steinStatus == 0 || seams_last_error() == nullptr)) {
      return fail("steinhardt guard");
    }
    return 0;
  };
  bad += probe_arm(arm77, "steinhardt guard");
  bad += arm77(false) != 0;
  }
  const std::string missingXyz = (std::filesystem::temp_directory_path() / "dseams-lang-missing.xyz").string();
    {
  auto arm78 = [&](bool force) -> int {
    if (force || (seams_frame_read(missingXyz.c_str(), 1, -1) != nullptr)) {
      return fail("missing xyz should fail");
    }
    return 0;
  };
  bad += probe_arm(arm78, "missing xyz should fail");
  bad += arm78(false) != 0;
  }
    {
  auto arm79 = [&](bool force) -> int {
    if (force || (seams_frame_chill(tip, 3.5, 99, "cutoff", labels.data(), counts) != 0)) {
      return fail("chill type filter");
    }
    return 0;
  };
  bad += probe_arm(arm79, "chill type filter");
  bad += arm79(false) != 0;
  }
    {
  auto arm80 = [&](bool force) -> int {
    if (force || (seams_frame_cages(built, 3.0, 99, 4, "nope", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) == 0)) {
      return fail("cages should reject the graph");
    }
    return 0;
  };
  bad += probe_arm(arm80, "cages should reject the graph");
  bad += arm80(false) != 0;
  }
    {
  auto arm81 = [&](bool force) -> int {
    if (force || (seams_frame_cages(built, 3.0, 99, 4, "cutoff", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) != 0)) {
      return fail("cages type filter");
    }
    return 0;
  };
  bad += probe_arm(arm81, "cages type filter");
  bad += arm81(false) != 0;
  }
  seams_prism_row oneRow{};
  int shortRows = 1;
    {
  auto arm82 = [&](bool force) -> int {
    if (force || (seams_frame_prisms(built, 3.0, 1, 6, nullptr, &oneRow, &shortRows) == 0)) {
      return fail("short prism buffer");
    }
    return 0;
  };
  bad += probe_arm(arm82, "short prism buffer");
  bad += arm82(false) != 0;
  }
    {
  auto arm83 = [&](bool force) -> int {
    if (force || (seams_frame_steinhardt(built, 3.0, 1, -1, qlEdge, qlbarEdge) != 0 && seams_last_error() == nullptr)) {
      return fail("steinhardt error");
    }
    return 0;
  };
  bad += probe_arm(arm83, "steinhardt error");
  bad += arm83(false) != 0;
  }
    {
  auto arm84 = [&](bool force) -> int {
    if (force || (seams_frame_rdf(built, -1.0, 4, 1, 1, radius, gram, histEdge, nullptr) != 0 && seams_last_error() == nullptr)) {
      return fail("rdf error");
    }
    return 0;
  };
  bad += probe_arm(arm84, "rdf error");
  bad += arm84(false) != 0;
  }
    {
  auto arm85 = [&](bool force) -> int {
    if (force || (seams_frame_rdf2d(built, -1.0, 0.1, nullptr, radius, gram) != 0 && seams_last_error() == nullptr)) {
      return fail("rdf2d error");
    }
    return 0;
  };
  bad += probe_arm(arm85, "rdf2d error");
  bad += arm85(false) != 0;
  }
    {
  auto arm86 = [&](bool force) -> int {
    if (force || (seams_frame_monolayer(built, -1.0, 1, 6, 1.0, nullptr, mono, 9) != 0 && seams_last_error() == nullptr)) {
      return fail("monolayer error");
    }
    return 0;
  };
  bad += probe_arm(arm86, "monolayer error");
  bad += arm86(false) != 0;
  }
  std::string longXyz(200000, 'q');
  longXyz += ".lammpstrj";
  if (seams_path_nframes(longXyz.c_str()) != -1) {
    return fail("long path frames");
  }
  const std::string oneAtom =
      (std::filesystem::temp_directory_path() / "dseams-lang-one.xyz").string();
  {
    std::ofstream out(oneAtom);
    out << "1\n\nO 0 0 0\n";
  }
  seams_frame *one = seams_frame_read(oneAtom.c_str(), 0, -1);
  if (one != nullptr) {
    double q = 0;
    double qb = 0;
    if (seams_frame_steinhardt(one, 3.0, 1, 4, &q, &qb) != 0) {
      return fail("one atom steinhardt");
    }
    if (seams_frame_steinhardt(one, 3.0, 1, 0, &q, &qb) != 0) {
      return fail("one atom steinhardt order");
    }
    seams_frame_free(one);
  }
  std::filesystem::remove(oneAtom);
  const auto blocker =
      std::filesystem::temp_directory_path() / "dseams-lang-notdir";
  {
    std::ofstream out(blocker);
    out << "x\n";
  }
  seams_prism_row prismRows[8]{};
  int prismN = 8;
  if (seams_frame_prisms(built, 3.0, 1, 6, blocker.string().c_str(),
                         prismRows, &prismN) != 0) {
    return fail("prism directory");
  }
  int monoCounts[8] = {};
  if (seams_frame_monolayer(built, 3.0, 1, 6, 1.0, blocker.string().c_str(),
                            monoCounts, 8) != 0) {
    return fail("monolayer directory");
  }
  std::filesystem::remove(blocker);

  seams_frame_free(fromXyz);
  seams_frame_free(empty);
  seams_frame_free(built);
  seams_frame_free(tilted);
  seams_frame_free(nullptr);

  seams_frame_free(molten);
  seams_frame_free(tip);
  seams_frame_free(ice);
  seams_frame_free(guest);
  if (bad != 0) {
    return 1;
  }
  fail("missing trajectory");
  std::cout << "seams_lang ok " << seams_lang_version() << "\n";
  return 0;
}
