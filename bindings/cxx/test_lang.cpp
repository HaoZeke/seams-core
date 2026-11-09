#include "seams_lang.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <cstring>
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

// A false condition prints the message and returns 1. The self-test
// below calls check(false) and reads that text back.
int check(bool ok, const char *message) {
  if (!ok) {
    return fail(message);
  }
  return 0;
}



int self_test_fail_arm() {
  struct rlimit nofile {};
  getrlimit(RLIMIT_NOFILE, &nofile);
  struct rlimit tight = nofile;
  tight.rlim_cur = 1;
  if (setrlimit(RLIMIT_NOFILE, &tight) == 0) {
    int blocked[2];
    if (pipe(blocked) != 0) {
      if (check(false, "fail arm pipe") == 0) {
        setrlimit(RLIMIT_NOFILE, &nofile);
        return 1;
      }
    } else {
      close(blocked[0]);
      close(blocked[1]);
    }
    setrlimit(RLIMIT_NOFILE, &nofile);
  }
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
  int rc = 0;
  rc = check(status == 1, "fail arm status");
  rc = rc != 0 ? rc : check(text.find("fail arm") != std::string::npos, "fail arm text");
  rc = rc != 0 ? rc : check(text.find("path is null") != std::string::npos, "fail arm error");
  return rc;
}

int count_of(const int counts[9], int label) { return counts[label]; }

} // namespace

int main(int argc, char **argv) {
  int rc = 0;
  const std::string root = argc > 1 ? argv[1] : ".";
  const std::string tip4p = root + "/input/traj/exampleTraj.lammpstrj";
  const std::string cubic = root + "/input/traj/mW_cubic.lammpstrj";
  const std::string si = root + "/input/traj/genice_sI.lammpstrj";
  rc = self_test_fail_arm();
  if (rc != 0) {
    return rc;
  }

  seams_frame *tip = seams_frame_read(tip4p.c_str(), 1, 2);
    {
  auto arm4 = [&]() -> int {
    return check(!(tip == nullptr), "read tip4p");
  };
  rc = rc != 0 ? rc : arm4();
  }
    {
  auto arm5 = [&]() -> int {
    return check(!(seams_frame_nop(tip) != 250), "tip4p nop");
  };
  rc = rc != 0 ? rc : arm5();
  }
  const double *xyz = seams_frame_positions(tip);
    {
  auto arm6 = [&]() -> int {
    return check(!(xyz == nullptr), "positions");
  };
  rc = rc != 0 ? rc : arm6();
  }

  std::vector<int> labels(250);
  int counts[9] = {};
    {
  auto arm7 = [&]() -> int {
    return check(!(seams_frame_chill_plus(tip, 3.5, 2, "mutual", labels.data(), counts) != 0), "chill+ mutual");
  };
  rc = rc != 0 ? rc : arm7();
  }
    {
  auto arm8 = [&]() -> int {
    return check(!(count_of(counts, 4) != 1 || count_of(counts, 5) != 7 || count_of(counts, 2) != 242), "mutual chill+ counts");
  };
  rc = rc != 0 ? rc : arm8();

  }

    {
  auto arm9 = [&]() -> int {
    return check(!(seams_frame_chill_plus(tip, 3.5, 2, "cutoff", labels.data(), counts) != 0), "chill+ cutoff");
  };
  rc = rc != 0 ? rc : arm9();
  }
    {
  auto arm10 = [&]() -> int {
    return check(!(count_of(counts, 5) != 12 || count_of(counts, 2) != 238 || count_of(counts, 4) != 0), "cutoff chill+ counts");
  };
  rc = rc != 0 ? rc : arm10();

  }
  int default_counts[9] = {};
    {
  auto arm11 = [&]() -> int {
    return check(!(seams_frame_chill_plus(tip, 3.5, 2, nullptr, labels.data(), default_counts) != 0), "chill+ default");
  };
  rc = rc != 0 ? rc : arm11();
  }
  for (int i = 0; i < 9; ++i) {
        {
    auto arm12 = [&]() -> int {
      return check(!(default_counts[i] != counts[i]), "chill+ default is not the cutoff list");
    };
    rc = rc != 0 ? rc : arm12();

    }
  }
    {
  auto arm13 = [&]() -> int {
    return check(seams_frame_chill_plus(tip, 3.5, 2, "knn-union", labels.data(), counts) != 0 &&
                     seams_last_error() != nullptr &&
                     std::string(seams_last_error()).find("unknown CHILL graph") != std::string::npos,
                 "chill+ accepted knn-union");
  };
  rc = rc != 0 ? rc : arm13();
  }

  int n_hex = -1, n_cubic = -1, n_water = -1, n_both = -1;
    {
  auto arm14 = [&]() -> int {
    return check(!(seams_frame_cages(tip, 3.5, 2, 4, "seeded", 0, nullptr, &n_hex, &n_cubic, &n_water, &n_both) != 0), "cages");
  };
  rc = rc != 0 ? rc : arm14();
  }
    {
  auto arm15 = [&]() -> int {
    return check(!(n_hex != 0 || n_cubic != 0 || n_water != 250), "cage counts");
  };
  rc = rc != 0 ? rc : arm15();

  }

  std::vector<double> ql(250), qlbar(250);
    {
  auto arm16 = [&]() -> int {
    return check(!(seams_frame_steinhardt(tip, 3.5, 2, 3, ql.data(), qlbar.data()) != 0), "steinhardt");
  };
  rc = rc != 0 ? rc : arm16();
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
  auto arm17 = [&]() -> int {
    return check(!(std::abs(mean - 0.389362) > 1e-4), "q3");
  };
  rc = rc != 0 ? rc : arm17();
  }

  const int bins = 4;
  std::vector<double> r(bins), g(bins);
  std::vector<int> hist(bins);
  double volume = 0.0;
    {
  auto arm18 = [&]() -> int {
    return check(!(seams_frame_rdf(tip, 6.0, bins, 2, 2, r.data(), g.data(), hist.data(), &volume) != 0), "rdf");
  };
  rc = rc != 0 ? rc : arm18();
  }
    {
  auto arm19 = [&]() -> int {
    return check(!(hist[1] != 390 || hist[2] != 705 || hist[3] != 531), "rdf counts");
  };
  rc = rc != 0 ? rc : arm19();

  }

  std::vector<int> rings(9, 0);
    {
  auto arm20 = [&]() -> int {
    return check(!(seams_frame_rings(tip, 3.5, 2, "mutual", 8, rings.data(), 9) != 0), "rings");
  };
  rc = rc != 0 ? rc : arm20();
  }
  std::cout << "tip4p mutual rings";
  for (int i = 3; i <= 8; ++i) {
    std::cout << " " << i << ":" << rings[i];
  }
  std::cout << "\n";

  seams_frame *ice = seams_frame_read(cubic.c_str(), 1, 1);
    {
  auto arm21 = [&]() -> int {
    return check(!(ice == nullptr || seams_frame_nop(ice) != 4096), "mW");
  };
  rc = rc != 0 ? rc : arm21();
  }
  std::vector<int> ice_labels(4096);
  int ice_counts[9] = {};
    {
  auto arm22 = [&]() -> int {
    return check(!(seams_frame_chill_plus(ice, 3.5, 1, "mutual", ice_labels.data(), ice_counts) != 0), "mW chill+");
  };
  rc = rc != 0 ? rc : arm22();
  }
    {
  auto arm23 = [&]() -> int {
    return check(!(ice_counts[0] != 4096), "mW cubic count");
  };
  rc = rc != 0 ? rc : arm23();

  }
    {
  auto arm24 = [&]() -> int {
    return check(!(seams_frame_cages(ice, 3.5, 1, 4, "seeded", 0, nullptr, &n_hex, &n_cubic, &n_water, &n_both) != 0), "mW cages");
  };
  rc = rc != 0 ? rc : arm24();
  }
    {
  auto arm25 = [&]() -> int {
    return check(!(n_hex != 0 || n_cubic != 4096 || n_water != 0), "mW cage counts");
  };
  rc = rc != 0 ? rc : arm25();

  }
    {
  auto arm26 = [&]() -> int {
    return check(!(seams_path_nframes(cubic.c_str()) < 2), "mW nframes");
  };
  rc = rc != 0 ? rc : arm26();

  }

  seams_frame *guest = seams_frame_read(si.c_str(), 1, 1);
    {
  auto arm27 = [&]() -> int {
    return check(!(guest == nullptr), "sI");
  };
  rc = rc != 0 ? rc : arm27();
  }
  int n_cages = 0;
  int n_atoms = 0;
    {
  auto arm28 = [&]() -> int {
    return check(!(seams_frame_signature(guest, 3.5, 1, "cutoff", "512", &n_cages, &n_atoms) != 0), "sI signature");
  };
  rc = rc != 0 ? rc : arm28();
  }
    {
  auto arm29 = [&]() -> int {
    return check(!(seams_frame_nop(guest) != 46 || n_cages != 2 || n_atoms != 40), "sI 512 census");
  };
  rc = rc != 0 ? rc : arm29();
  }
    {
  auto arm30 = [&]() -> int {
    return check(!(seams_frame_signature_edge(guest, 3.5, 1, "cutoff", "512", "none", &n_cages, &n_atoms) != 0 || n_cages != 2 || n_atoms != 40), "sI 512 none");
  };
  rc = rc != 0 ? rc : arm30();
  }
    {
  auto arm31 = [&]() -> int {
    return check(!(seams_frame_signature_edge(guest, 3.5, 1, "cutoff", "512", "1,2", &n_cages, &n_atoms) != 0 || n_cages != 0), "sI 512 tight window");
  };
  rc = rc != 0 ? rc : arm31();
  }
    {
  auto arm32 = [&]() -> int {
    return check(!(seams_frame_signature(guest, 3.5, 1, "cutoff", "51268", &n_cages, &n_atoms) != 0), "sI 51268");
  };
  rc = rc != 0 ? rc : arm32();
  }
    {
  auto arm33 = [&]() -> int {
    return check(!(n_cages != 0), "sI 51268 census");
  };
  rc = rc != 0 ? rc : arm33();
  }

  const std::string nucleation = root + "/input/traj/nucleation_f1.lammpstrj";
  seams_frame *molten = seams_frame_read(nucleation.c_str(), 1, 1);
    {
  auto arm34 = [&]() -> int {
    return check(!(molten == nullptr), "nucleation frame");
  };
  rc = rc != 0 ? rc : arm34();
  }
  const int nMolten = seams_frame_nop(molten);
  std::vector<int> hostLabels(static_cast<std::size_t>(nMolten));
  std::vector<int> offLabels(static_cast<std::size_t>(nMolten));
  int hHex = 0, hCub = 0, hWat = 0, hBoth = 0;
  int oHex = 0, oCub = 0, oWat = 0, oBoth = 0;
  unsetenv("SEAMS_OFFLOAD");
    {
  auto arm35 = [&]() -> int {
    return check(!(seams_frame_cages(molten, 3.5, 1, 4, "seeded", 0, hostLabels.data(), &hHex, &hCub, &hWat, &hBoth) != 0), "seeded host");
  };
  rc = rc != 0 ? rc : arm35();
  }
  setenv("SEAMS_OFFLOAD", "1", 1);
    {
  auto arm36 = [&]() -> int {
    return check(!(seams_frame_cages(molten, 3.5, 1, 4, "seeded", 0, offLabels.data(), &oHex, &oCub, &oWat, &oBoth) != 0), "seeded offload");
  };
  rc = rc != 0 ? rc : arm36();
  }
  unsetenv("SEAMS_OFFLOAD");
    {
  auto arm37 = [&]() -> int {
    return check(!(hHex != oHex || hCub != oCub || hWat != oWat || hBoth != oBoth || hostLabels != offLabels), "seeded offload disagrees with the host rule");
  };
  rc = rc != 0 ? rc : arm37();
  }
    {
  auto arm38 = [&]() -> int {
    return check(!(hHex + hCub + hBoth + hWat != nMolten), "disjoint cage counters");
  };
  rc = rc != 0 ? rc : arm38();
  }
  int codes[4] = {};
  for (int code : hostLabels) {
        {
    auto arm39 = [&]() -> int {
      return check(!(code < 0 || code > 3), "cage code");
    };
    rc = rc != 0 ? rc : arm39();
    }
    ++codes[code];
  }
  // 0 water, 1 hexagonal only, 2 cubic only, 3 both
    {
  auto arm40 = [&]() -> int {
    return check(!(codes[1] != hHex || codes[2] != hCub || codes[3] != hBoth || codes[0] != hWat || hBoth <= 0), "stacking-fault counters");
  };
  rc = rc != 0 ? rc : arm40();

  }
  const std::string cli = argc > 2 ? std::string(argv[2]) : std::string("bbdir-cov/seams");
  if (!cli.empty()) {
    const std::string cmd = cli + " cages " + nucleation +
                            " --graph seeded --cutoff 3.5";
    FILE *pipe = popen(cmd.c_str(), "r");
        {
    auto arm41 = [&]() -> int {
      return check(!(pipe == nullptr), "cli pipe");
    };
    rc = rc != 0 ? rc : arm41();
    }
    std::string line;
    char buf[512];
    while (fgets(buf, sizeof buf, pipe) != nullptr) {
      line += buf;
    }
    const int cli_rc = pclose(pipe);
        {
    auto arm42 = [&]() -> int {
      return check(!(cli_rc != 0), "cli cages");
    };
    rc = rc != 0 ? rc : arm42();
    }
    int cliNop = 0, cliHex = 0, cliCub = 0, cliBoth = 0, cliWat = 0;
    std::istringstream in(line);
    std::string tok;
    in >> tok >> cliNop >> tok >> tok >> tok >> cliHex >> tok >> cliCub >>
        tok >> cliBoth >> tok >> cliWat;
        {
    auto arm43 = [&]() -> int {
      return check(!(cliHex != hHex || cliCub != hCub || cliBoth != hBoth || cliWat != hWat), "abi disagrees with the cli");
    };
    rc = rc != 0 ? rc : arm43();
    }
  }
  {
 auto arm44 = [&]() -> int {
   return check(!(seams_last_error() == nullptr), "error pointer");
 };
 rc = rc != 0 ? rc : arm44();
 }
    {
  auto arm45 = [&]() -> int {
    return check(!(seams_path_nframes(nullptr) != -1 || seams_frame_read(nullptr, 1, -1) != nullptr), "null path");
  };
  rc = rc != 0 ? rc : arm45();
  }
  double box9[9] = {};
  seams_frame_box(nullptr, box9);
  seams_frame_box(tip, nullptr);
  seams_frame_box(tip, box9);
    {
  auto arm46 = [&]() -> int {
    return check(!(seams_frame_index(nullptr) != 0 || seams_frame_nop(nullptr) != 0), "null frame queries");
  };
  rc = rc != 0 ? rc : arm46();
  }
    {
  auto arm47 = [&]() -> int {
    return check(!(seams_frame_index(tip) < 1 || seams_frame_types(tip) == nullptr || seams_frame_ids(tip) == nullptr || seams_frame_types(nullptr) != nullptr || seams_frame_ids(nullptr) != nullptr || seams_frame_positions(nullptr) != nullptr), "frame columns");
  };
  rc = rc != 0 ? rc : arm47();
  }
    {
  auto arm48 = [&]() -> int {
    return check(!(seams_frame_from_arrays(nullptr, 4, box9, nullptr, nullptr) != nullptr), "null arrays");
  };
  rc = rc != 0 ? rc : arm48();
  }
  const double builtXyz[12] = {0, 0, 0, 1.2, 0, 0, 0, 1.2, 0, 0.4, 0.4, 1.2};
  const double ortho[9] = {0, 0, 0, 8, 8, 8, 0, 0, 0};
  const double sheared[9] = {0, 0, 0, 8, 8, 8, 1.5, 0, 0};
  const int types4[4] = {1, 1, 1, 1};
  seams_frame *built = seams_frame_from_arrays(builtXyz, 4, ortho, types4, nullptr);
  seams_frame *tilted = seams_frame_from_arrays(builtXyz, 4, sheared, nullptr, nullptr);
    {
  auto arm49 = [&]() -> int {
    return check(!(built == nullptr || tilted == nullptr), "from arrays");
  };
  rc = rc != 0 ? rc : arm49();
  }
  int labels4[4] = {};
  int counts9[9] = {};
    {
  auto arm50 = [&]() -> int {
    return check(!(seams_frame_chill(built, 3.5, 1, "cutoff", labels4, counts9) != 0 || seams_frame_chill(nullptr, 3.5, 1, "mutual", labels4, counts9) == 0 || seams_frame_chill(built, 3.5, 0, "nope", labels4, counts9) == 0), "chill");
  };
  rc = rc != 0 ? rc : arm50();
  }
  int ringCounts[9] = {};
    {
  auto arm51 = [&]() -> int {
    return check(!(seams_frame_rings(built, 3.0, 1, "knn-union", 6, ringCounts, 9) != 0 || seams_frame_rings(built, 3.0, 1, "cutoff", 6, ringCounts, 9) != 0 || seams_frame_rings(built, 3.0, 1, "knn", 6, ringCounts, 9) != 0 || seams_frame_rings(built, 3.0, 0, nullptr, 6, ringCounts, 9) != 0 || seams_frame_rings(nullptr, 3.0, 1, "mutual", 6, ringCounts, 9) == 0 || seams_frame_rings(built, 3.0, 1, "seeded", 6, ringCounts, 9) == 0 || seams_frame_rings(built, 3.0, 1, "nope", 6, ringCounts, 9) == 0), "rings edges");
  };
  rc = rc != 0 ? rc : arm51();
  }
    {
  auto arm52 = [&]() -> int {
    return check(!(seams_frame_rdf_nbins(6.0, 0) < 1 || seams_frame_rdf2d_nbin(3.0, 0.0) < 1 || seams_frame_rdf2d_nbin(0.01, 1.0) < 1), "bin counts");
  };
  rc = rc != 0 ? rc : arm52();
  }
  // nbins 0 asks for rmax/0.1 bins. 4.0 / 0.1 is 40.
  double radius[64] = {};
  double gram[64] = {};
  int histEdge[64] = {};
  double volumeEdge = 0.0;
    {
  auto arm53 = [&]() -> int {
    return check(!(seams_frame_rdf(built, 4.0, 0, 1, 0, radius, gram, histEdge, &volumeEdge) != 0 || seams_frame_rdf(nullptr, 4.0, 4, 1, 1, radius, gram, histEdge, nullptr) == 0), "rdf edges");
  };
  rc = rc != 0 ? rc : arm53();
  }
  const std::string outdir = "var/run/dseams-lang-out";
    {
  auto arm54 = [&]() -> int {
    return check(!(seams_frame_rdf2d(built, 3.0, 0.5, nullptr, radius, gram) != 0 || seams_frame_rdf2d(built, 3.0, 0.5, outdir.c_str(), radius, gram) != 0 || seams_frame_rdf2d(nullptr, 3.0, 0.5, nullptr, radius, gram) == 0), "rdf2d");
  };
  rc = rc != 0 ? rc : arm54();
  }
  int nrows = 0;
    {
  auto arm55 = [&]() -> int {
    return check(!(seams_frame_prisms(built, 3.0, 1, 6, nullptr, nullptr, &nrows) != 0 || nrows < 1), "prism count");
  };
  rc = rc != 0 ? rc : arm55();
  }
  std::vector<seams_prism_row> prows(static_cast<std::size_t>(nrows));
    {
  auto arm56 = [&]() -> int {
    return check(!(seams_frame_prisms(nullptr, 3.0, 1, 2, nullptr, nullptr, nullptr) == 0 || seams_frame_prisms(built, 3.0, 1, 6, outdir.c_str(), prows.data(), &nrows) != 0), "prisms");
  };
  rc = rc != 0 ? rc : arm56();
  }
  int mono[9] = {};
    {
  auto arm57 = [&]() -> int {
    return check(!(seams_frame_monolayer(built, 3.0, 1, 6, 10.0, nullptr, mono, 9) != 0 || seams_frame_monolayer(built, 3.0, 1, 6, 10.0, outdir.c_str(), mono, 9) != 0 || seams_frame_monolayer(nullptr, 3.0, 1, 6, 10.0, nullptr, mono, 9) == 0), "monolayer");
  };
  rc = rc != 0 ? rc : arm57();
  }
  int nHex = 0, nCub = 0, nWat = 0, nBoth = 0;
  int perAtom[4] = {};
    {
  auto arm58 = [&]() -> int {
    return check(!(seams_frame_cages(built, 3.0, 1, 4, "cutoff", 1, perAtom, &nHex, &nCub, &nWat, &nBoth) != 0 || seams_frame_cages(nullptr, 3.0, 1, 4, "seeded", 0, nullptr, &nHex, &nCub, &nWat, &nBoth) == 0), "cage edges");
  };
  rc = rc != 0 ? rc : arm58();
  }
  double qlEdge[4] = {};
  double qlbarEdge[4] = {};
    {
  auto arm59 = [&]() -> int {
    return check(!(seams_frame_steinhardt(nullptr, 3.0, 1, 4, qlEdge, qlbarEdge) == 0 || seams_frame_steinhardt(built, 3.0, 1, 4, nullptr, nullptr) != 0), "steinhardt edges");
  };
  rc = rc != 0 ? rc : arm59();
  }
  int nCages = 0, nAtoms = 0;
    {
  auto arm60 = [&]() -> int {
    return check(seams_frame_signature(nullptr, 3.0, 1, "cutoff", "512", &nCages, &nAtoms) != 0 &&
                     seams_frame_signature(built, 3.0, 1, "nope", "512", &nCages, &nAtoms) != 0 &&
                     seams_last_error() != nullptr &&
                     std::string(seams_last_error()).find("unknown bond graph") != std::string::npos,
                 "signature edges");
  };
  rc = rc != 0 ? rc : arm60();
  }
  double tiltedBox[9] = {};
  seams_frame_box(tilted, tiltedBox);
    {
  auto arm61 = [&]() -> int {
    return check(!(tiltedBox[6] == 0.0), "shear was dropped");
  };
  rc = rc != 0 ? rc : arm61();
  }
    {
  auto arm62 = [&]() -> int {
    return check(!(seams_path_nframes("nodot") < 0 && seams_last_error() == nullptr), "extension");
  };
  rc = rc != 0 ? rc : arm62();
  }
  int guestHex = 0, guestCub = 0, guestWat = 0, guestBoth = 0;
    {
  auto arm63 = [&]() -> int {
    return check(!(seams_frame_cages(guest, 3.5, 1, 0, "cutoff", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) != 0), "guest cutoff cages");
  };
  rc = rc != 0 ? rc : arm63();
  }
  setenv("SEAMS_OFFLOAD", "1", 1);
    {
  auto arm64 = [&]() -> int {
    return check(!(seams_frame_cages(guest, 3.5, 1, 4, "mutual", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) != 0), "guest offload cages");
  };
  rc = rc != 0 ? rc : arm64();
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
  auto arm65 = [&]() -> int {
    return check(!(empty == nullptr || seams_frame_cages(empty, 3.0, 0, 4, "cutoff", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) != 0), "empty cages");
  };
  rc = rc != 0 ? rc : arm65();
  }
  const std::string tinyXyz = (std::filesystem::temp_directory_path() / "dseams-lang-tiny.xyz").string();
  {
    std::ofstream out(tinyXyz);
    out << "1\ncomment\nO 0 0 0\n";
  }
  seams_frame *fromXyz = seams_frame_read(tinyXyz.c_str(), 1, -1);
    {
  auto arm66 = [&]() -> int {
    return check(!(fromXyz == nullptr), "xyz read");
  };
  rc = rc != 0 ? rc : arm66();
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
  auto arm67 = [&]() -> int {
    return check(!(seams_path_nframes(twoXyz.c_str()) != 2), "xyz frame count");
  };
  rc = rc != 0 ? rc : arm67();

  }
    {
  auto arm68 = [&]() -> int {
    return check(!(seams_path_nframes(oneGro.c_str()) != 1), "gro frame count");
  };
  rc = rc != 0 ? rc : arm68();

  }
    {
  auto arm69 = [&]() -> int {
    return check(!(seams_path_nframes(onePdb.c_str()) != 1), "pdb frame count");
  };
  rc = rc != 0 ? rc : arm69();

  }
  const std::string missingXyzCount = (nroot / "dseams-nframes-missing.xyz").string();
  const std::string badXyz = (nroot / "dseams-nframes-bad.xyz").string();
  {
    std::ofstream out(badXyz);
    out << "this is not an xyz file\n";
  }
    {
  auto arm70 = [&]() -> int {
    return check(!(seams_path_nframes(missingXyzCount.c_str()) != -1 || seams_last_error() == nullptr), "missing xyz count");
  };
  rc = rc != 0 ? rc : arm70();

  }
    {
  auto arm71 = [&]() -> int {
    return check(!(seams_path_nframes(badXyz.c_str()) != -1 || seams_last_error() == nullptr), "bad xyz count");
  };
  rc = rc != 0 ? rc : arm71();

  }
    {
  auto arm72 = [&]() -> int {
    return check(!(seams_path_nframes("missing-file.lammpstrj") != -1 || seams_last_error() == nullptr), "missing dump count");
  };
  rc = rc != 0 ? rc : arm72();
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
    setrlimit(RLIMIT_AS, &lim);
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
  auto arm73 = [&]() -> int {
    return check(!(over != nullptr || seams_last_error() == nullptr), "array guard");
  };
  rc = rc != 0 ? rc : arm73();
  }
  seams_frame_free(over);
  std::vector<double> rr(static_cast<std::size_t>(many));
  std::vector<double> gg(static_cast<std::size_t>(many));
  std::vector<int> cc(static_cast<std::size_t>(many));
  cap = tighten();
  const int rdfStatus = seams_frame_rdf(built, 10.0, many, 1, 1, rr.data(), gg.data(), cc.data(), nullptr);
  loosen(cap);
    {
  auto arm74 = [&]() -> int {
    return check(!(rdfStatus == 0 || seams_last_error() == nullptr), "rdf guard");
  };
  rc = rc != 0 ? rc : arm74();
  }
  cap = tighten();
  const int rdf2Status = seams_frame_rdf2d(built, 1.0e8, 0.05, nullptr, rr.data(), gg.data());
  loosen(cap);
    {
  auto arm75 = [&]() -> int {
    return check(!(rdf2Status == 0 || seams_last_error() == nullptr), "rdf2d guard");
  };
  rc = rc != 0 ? rc : arm75();
  }
  seams_frame *crowd = seams_frame_from_arrays(wide.data(), many, wideBox, nullptr, nullptr);
    {
  auto arm76 = [&]() -> int {
    return check(!(crowd == nullptr), "crowd frame");
  };
  rc = rc != 0 ? rc : arm76();
  }
  cap = tighten();
  const int steinStatus = seams_frame_steinhardt(crowd, 3.5, 1, 4, rr.data(), gg.data());
  loosen(cap);
  seams_frame_free(crowd);
    {
  auto arm77 = [&]() -> int {
    return check(!(steinStatus == 0 || seams_last_error() == nullptr), "steinhardt guard");
  };
  rc = rc != 0 ? rc : arm77();
  }
  const std::string missingXyz = (std::filesystem::temp_directory_path() / "dseams-lang-missing.xyz").string();
    {
  auto arm78 = [&]() -> int {
    return check(!(seams_frame_read(missingXyz.c_str(), 1, -1) != nullptr), "missing xyz should fail");
  };
  rc = rc != 0 ? rc : arm78();

  }
    {
  auto arm79 = [&]() -> int {
    return check(!(seams_frame_chill(tip, 3.5, 99, "cutoff", labels.data(), counts) != 0), "chill type filter");
  };
  rc = rc != 0 ? rc : arm79();
  }
    {
  auto arm80 = [&]() -> int {
    return check(seams_frame_cages(built, 3.0, 99, 4, "nope", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) != 0 &&
                     seams_last_error() != nullptr &&
                     std::string(seams_last_error()).find("unknown") != std::string::npos,
                 "cages should reject the graph");
  };
  rc = rc != 0 ? rc : arm80();
  }
    {
  auto arm81 = [&]() -> int {
    return check(!(seams_frame_cages(built, 3.0, 99, 4, "cutoff", 0, nullptr, &guestHex, &guestCub, &guestWat, &guestBoth) != 0), "cages type filter");
  };
  rc = rc != 0 ? rc : arm81();
  }
  seams_prism_row oneRow{};
  int shortRows = 1;
    {
  auto arm82 = [&]() -> int {
    return check(seams_frame_prisms(built, 3.0, 1, 6, nullptr, &oneRow, &shortRows) != 0 &&
                     seams_last_error() != nullptr &&
                     std::string(seams_last_error()).find("prism row buffer is short") != std::string::npos,
                 "short prism buffer");
  };
  rc = rc != 0 ? rc : arm82();
  }
    {
  auto arm83 = [&]() -> int {
    const int steinEdge = seams_frame_steinhardt(built, 3.0, 1, -1, qlEdge, qlbarEdge);
    return check(!(steinEdge != 0 && (seams_last_error() == nullptr || seams_last_error()[0] == '\0')),
                 "steinhardt error");
  };
  rc = rc != 0 ? rc : arm83();
  }
    {
  auto arm84 = [&]() -> int {
    const int rdfEdge = seams_frame_rdf(built, -1.0, 4, 1, 1, radius, gram, histEdge, nullptr);
    return check(!(rdfEdge != 0 && (seams_last_error() == nullptr || seams_last_error()[0] == '\0')),
                 "rdf error");
  };
  rc = rc != 0 ? rc : arm84();
  }
    {
  auto arm85 = [&]() -> int {
    const int rdf2Edge = seams_frame_rdf2d(built, -1.0, 0.1, nullptr, radius, gram);
    return check(!(rdf2Edge != 0 && (seams_last_error() == nullptr || seams_last_error()[0] == '\0')),
                 "rdf2d error");
  };
  rc = rc != 0 ? rc : arm85();
  }
    {
  auto arm86 = [&]() -> int {
    const int monoEdge = seams_frame_monolayer(built, -1.0, 1, 6, 1.0, nullptr, mono, 9);
    return check(!(monoEdge != 0 && (seams_last_error() == nullptr || seams_last_error()[0] == '\0')),
                 "monolayer error");
  };
  rc = rc != 0 ? rc : arm86();
  }
  std::string longXyz(200000, 'q');
  longXyz += ".lammpstrj";
  rc = rc != 0 ? rc : check(seams_path_nframes(longXyz.c_str()) == -1, "long path frames");
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
    rc = rc != 0 ? rc : check(seams_frame_steinhardt(one, 3.0, 1, 4, &q, &qb) == 0, "one atom steinhardt");
    rc = rc != 0 ? rc : check(seams_frame_steinhardt(one, 3.0, 1, 0, &q, &qb) == 0, "one atom steinhardt order");
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
  rc = rc != 0 ? rc : check(seams_frame_prisms(built, 3.0, 1, 6, blocker.string().c_str(),
                                              prismRows, &prismN) == 0,
                            "prism directory");
  int monoCounts[8] = {};
  rc = rc != 0 ? rc : check(seams_frame_monolayer(built, 3.0, 1, 6, 1.0, blocker.string().c_str(),
                                                  monoCounts, 8) == 0,
                            "monolayer directory");
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
  if (rc != 0) {
    return rc;
  }
  std::cout << "seams_lang ok " << seams_lang_version() << "\n";
  return rc;
}
