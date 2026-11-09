#include "seams_lang.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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

int count_of(const int counts[9], int label) { return counts[label]; }

} // namespace

int main(int argc, char **argv) {
  const std::string root = argc > 1 ? argv[1] : ".";
  const std::string tip4p = root + "/input/traj/exampleTraj.lammpstrj";
  const std::string cubic = root + "/input/traj/mW_cubic.lammpstrj";
  const std::string si = root + "/input/traj/genice_sI.lammpstrj";

  if (std::string(seams_lang_version()).empty()) {
    return fail("empty version");
  }

  seams_frame *tip = seams_frame_read(tip4p.c_str(), 1, 2);
  if (tip == nullptr) {
    return fail("read tip4p");
  }
  if (seams_frame_nop(tip) != 250) {
    return fail("tip4p nop");
  }
  const double *xyz = seams_frame_positions(tip);
  if (xyz == nullptr) {
    return fail("positions");
  }

  std::vector<int> labels(250);
  int counts[9] = {};
  if (seams_frame_chill_plus(tip, 3.5, 2, "mutual", labels.data(), counts) !=
      0) {
    return fail("chill+ mutual");
  }
  if (count_of(counts, 4) != 1 || count_of(counts, 5) != 7 ||
      count_of(counts, 2) != 242) {
    std::cerr << "mutual chill+ " << counts[4] << " " << counts[5] << " "
              << counts[2] << "\n";
    return fail("mutual chill+ counts");
  }

  if (seams_frame_chill_plus(tip, 3.5, 2, "cutoff", labels.data(), counts) !=
      0) {
    return fail("chill+ cutoff");
  }
  if (count_of(counts, 5) != 12 || count_of(counts, 2) != 238 ||
      count_of(counts, 4) != 0) {
    std::cerr << "cutoff chill+ clathrate " << counts[4] << " inter "
              << counts[5] << " water " << counts[2] << "\n";
    return fail("cutoff chill+ counts");
  }
  int default_counts[9] = {};
  if (seams_frame_chill_plus(tip, 3.5, 2, nullptr, labels.data(),
                             default_counts) != 0) {
    return fail("chill+ default");
  }
  for (int i = 0; i < 9; ++i) {
    if (default_counts[i] != counts[i]) {
      return fail("chill+ default is not the cutoff list");
    }
  }
  if (seams_frame_chill_plus(tip, 3.5, 2, "knn-union", labels.data(),
                             counts) == 0) {
    return fail("chill+ accepted knn-union");
  }

  int n_hex = -1, n_cubic = -1, n_water = -1, n_both = -1;
  if (seams_frame_cages(tip, 3.5, 2, 4, "seeded", 0, nullptr, &n_hex, &n_cubic,
                        &n_water, &n_both) != 0) {
    return fail("cages");
  }
  if (n_hex != 0 || n_cubic != 0 || n_water != 250) {
    std::cerr << "cages " << n_hex << " " << n_cubic << " " << n_water << "\n";
    return fail("cage counts");
  }

  std::vector<double> ql(250), qlbar(250);
  if (seams_frame_steinhardt(tip, 3.5, 2, 3, ql.data(), qlbar.data()) != 0) {
    return fail("steinhardt");
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
  if (std::abs(mean - 0.389362) > 1e-4) {
    std::cerr << "q3 mean " << mean << "\n";
    return fail("q3");
  }

  const int bins = 4;
  std::vector<double> r(bins), g(bins);
  std::vector<int> hist(bins);
  double volume = 0.0;
  if (seams_frame_rdf(tip, 6.0, bins, 2, 2, r.data(), g.data(), hist.data(),
                      &volume) != 0) {
    return fail("rdf");
  }
  if (hist[1] != 390 || hist[2] != 705 || hist[3] != 531) {
    std::cerr << "rdf counts " << hist[0] << " " << hist[1] << " " << hist[2]
              << " " << hist[3] << "\n";
    return fail("rdf counts");
  }

  std::vector<int> rings(9, 0);
  if (seams_frame_rings(tip, 3.5, 2, "mutual", 8, rings.data(), 9) != 0) {
    return fail("rings");
  }
  std::cout << "tip4p mutual rings";
  for (int i = 3; i <= 8; ++i) {
    std::cout << " " << i << ":" << rings[i];
  }
  std::cout << "\n";

  seams_frame *ice = seams_frame_read(cubic.c_str(), 1, 1);
  if (ice == nullptr || seams_frame_nop(ice) != 4096) {
    return fail("mW");
  }
  std::vector<int> ice_labels(4096);
  int ice_counts[9] = {};
  if (seams_frame_chill_plus(ice, 3.5, 1, "mutual", ice_labels.data(),
                             ice_counts) != 0) {
    return fail("mW chill+");
  }
  if (ice_counts[0] != 4096) {
    std::cerr << "mW cubic " << ice_counts[0] << "\n";
    return fail("mW cubic count");
  }
  if (seams_frame_cages(ice, 3.5, 1, 4, "seeded", 0, nullptr, &n_hex, &n_cubic,
                        &n_water, &n_both) != 0) {
    return fail("mW cages");
  }
  if (n_hex != 0 || n_cubic != 4096 || n_water != 0) {
    std::cerr << "mW cages " << n_hex << " " << n_cubic << " " << n_water
              << "\n";
    return fail("mW cage counts");
  }
  if (seams_path_nframes(cubic.c_str()) < 2) {
    return fail("mW nframes");
  }

  seams_frame *guest = seams_frame_read(si.c_str(), 1, 1);
  if (guest == nullptr) {
    return fail("sI");
  }
  int n_cages = 0;
  int n_atoms = 0;
  if (seams_frame_signature(guest, 3.5, 1, "cutoff", "512", &n_cages,
                            &n_atoms) != 0) {
    return fail("sI signature");
  }
  if (seams_frame_nop(guest) != 46 || n_cages != 2 || n_atoms != 40) {
    std::cerr << "sI " << seams_frame_nop(guest) << " cages " << n_cages
              << " atoms " << n_atoms << "\n";
    return fail("sI 512 census");
  }
  if (seams_frame_signature_edge(guest, 3.5, 1, "cutoff", "512", "none",
                                 &n_cages, &n_atoms) != 0 ||
      n_cages != 2 || n_atoms != 40) {
    std::cerr << "sI none cages " << n_cages << " atoms " << n_atoms << "\n";
    return fail("sI 512 none");
  }
  if (seams_frame_signature_edge(guest, 3.5, 1, "cutoff", "512", "1,2",
                                 &n_cages, &n_atoms) != 0 ||
      n_cages != 0) {
    std::cerr << "sI tight cages " << n_cages << "\n";
    return fail("sI 512 tight window");
  }
  if (seams_frame_signature(guest, 3.5, 1, "cutoff", "51268", &n_cages,
                            &n_atoms) != 0) {
    return fail("sI 51268");
  }
  if (n_cages != 0) {
    std::cerr << "sI 51268 cages " << n_cages << "\n";
    return fail("sI 51268 census");
  }

  const std::string nucleation = root + "/input/traj/nucleation_f1.lammpstrj";
  seams_frame *molten = seams_frame_read(nucleation.c_str(), 1, 1);
  if (molten == nullptr) {
    return fail("nucleation frame");
  }
  const int nMolten = seams_frame_nop(molten);
  std::vector<int> hostLabels(static_cast<std::size_t>(nMolten));
  std::vector<int> offLabels(static_cast<std::size_t>(nMolten));
  int hHex = 0, hCub = 0, hWat = 0, hBoth = 0;
  int oHex = 0, oCub = 0, oWat = 0, oBoth = 0;
  unsetenv("SEAMS_OFFLOAD");
  if (seams_frame_cages(molten, 3.5, 1, 4, "seeded", 0, hostLabels.data(),
                        &hHex, &hCub, &hWat, &hBoth) != 0) {
    return fail("seeded host");
  }
  setenv("SEAMS_OFFLOAD", "1", 1);
  if (seams_frame_cages(molten, 3.5, 1, 4, "seeded", 0, offLabels.data(),
                        &oHex, &oCub, &oWat, &oBoth) != 0) {
    return fail("seeded offload");
  }
  unsetenv("SEAMS_OFFLOAD");
  if (hHex != oHex || hCub != oCub || hWat != oWat || hBoth != oBoth ||
      hostLabels != offLabels) {
    std::cerr << "offload " << oHex << " " << oCub << " " << oBoth << " "
              << oWat << " host " << hHex << " " << hCub << " " << hBoth << " "
              << hWat << "\n";
    return fail("seeded offload disagrees with the host rule");
  }
  if (hHex + hCub + hBoth + hWat != nMolten) {
    std::cerr << "counters " << hHex << " " << hCub << " " << hBoth << " "
              << hWat << " nop " << nMolten << "\n";
    return fail("disjoint cage counters");
  }
  int codes[4] = {};
  for (int code : hostLabels) {
    if (code < 0 || code > 3) {
      return fail("cage code");
    }
    ++codes[code];
  }
  // 0 water, 1 hexagonal only, 2 cubic only, 3 both
  if (codes[1] != hHex || codes[2] != hCub || codes[3] != hBoth ||
      codes[0] != hWat || hBoth <= 0) {
    std::cerr << "codes " << codes[1] << " " << codes[2] << " " << codes[3]
              << " both " << hBoth << "\n";
    return fail("stacking-fault counters");
  }
  if (argc > 2) {
    const std::string cmd = std::string(argv[2]) + " cages " + nucleation +
                            " --graph seeded --cutoff 3.5";
    FILE *pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr) {
      return fail("cli pipe");
    }
    std::string line;
    char buf[512];
    while (fgets(buf, sizeof buf, pipe) != nullptr) {
      line += buf;
    }
    const int rc = pclose(pipe);
    if (rc != 0) {
      std::cerr << line;
      return fail("cli cages");
    }
    int cliNop = 0, cliHex = 0, cliCub = 0, cliBoth = 0, cliWat = 0;
    std::istringstream in(line);
    std::string tok;
    in >> tok >> cliNop >> tok >> tok >> tok >> cliHex >> tok >> cliCub >>
        tok >> cliBoth >> tok >> cliWat;
    if (cliHex != hHex || cliCub != hCub || cliBoth != hBoth ||
        cliWat != hWat) {
      std::cerr << "cli " << line << "abi " << hHex << " " << hCub << " "
                << hBoth << " " << hWat << "\n";
      return fail("abi disagrees with the cli");
    }
  }

  if (seams_last_error() == nullptr) {
    return fail("error pointer");
  }
  if (seams_path_nframes(nullptr) != -1 || seams_frame_read(nullptr, 1, -1) != nullptr) {
    return fail("null path");
  }
  double box9[9] = {};
  seams_frame_box(nullptr, box9);
  seams_frame_box(tip, nullptr);
  seams_frame_box(tip, box9);
  if (seams_frame_index(nullptr) != 0 || seams_frame_nop(nullptr) != 0) {
    return fail("null frame queries");
  }
  if (seams_frame_index(tip) < 1 || seams_frame_types(tip) == nullptr ||
      seams_frame_ids(tip) == nullptr || seams_frame_types(nullptr) != nullptr ||
      seams_frame_ids(nullptr) != nullptr || seams_frame_positions(nullptr) != nullptr) {
    return fail("frame columns");
  }
  if (seams_frame_from_arrays(nullptr, 4, box9, nullptr, nullptr) != nullptr) {
    return fail("null arrays");
  }
  const double builtXyz[12] = {0, 0, 0, 1.2, 0, 0, 0, 1.2, 0, 0.4, 0.4, 1.2};
  const double ortho[9] = {0, 0, 0, 8, 8, 8, 0, 0, 0};
  const double sheared[9] = {0, 0, 0, 8, 8, 8, 1.5, 0, 0};
  const int types4[4] = {1, 1, 1, 1};
  seams_frame *built = seams_frame_from_arrays(builtXyz, 4, ortho, types4, nullptr);
  seams_frame *tilted = seams_frame_from_arrays(builtXyz, 4, sheared, nullptr, nullptr);
  if (built == nullptr || tilted == nullptr) {
    return fail("from arrays");
  }
  int labels4[4] = {};
  int counts9[9] = {};
  if (seams_frame_chill(built, 3.5, 1, "cutoff", labels4, counts9) != 0 ||
      seams_frame_chill(nullptr, 3.5, 1, "mutual", labels4, counts9) == 0 ||
      seams_frame_chill(built, 3.5, 0, "nope", labels4, counts9) == 0) {
    return fail("chill");
  }
  int ringCounts[9] = {};
  if (seams_frame_rings(built, 3.0, 1, "knn-union", 6, ringCounts, 9) != 0 ||
      seams_frame_rings(built, 3.0, 1, "cutoff", 6, ringCounts, 9) != 0 ||
      seams_frame_rings(built, 3.0, 1, "knn", 6, ringCounts, 9) != 0 ||
      seams_frame_rings(built, 3.0, 0, nullptr, 6, ringCounts, 9) != 0 ||
      seams_frame_rings(nullptr, 3.0, 1, "mutual", 6, ringCounts, 9) == 0 ||
      seams_frame_rings(built, 3.0, 1, "seeded", 6, ringCounts, 9) == 0 ||
      seams_frame_rings(built, 3.0, 1, "nope", 6, ringCounts, 9) == 0) {
    return fail("rings edges");
  }
  if (seams_frame_rdf_nbins(6.0, 0) < 1 || seams_frame_rdf2d_nbin(3.0, 0.0) < 1 ||
      seams_frame_rdf2d_nbin(0.01, 1.0) < 1) {
    return fail("bin counts");
  }
  double radius[8] = {};
  double gram[8] = {};
  int histEdge[8] = {};
  double volumeEdge = 0.0;
  if (seams_frame_rdf(built, 4.0, 0, 1, 0, radius, gram, histEdge, &volumeEdge) != 0 ||
      seams_frame_rdf(nullptr, 4.0, 4, 1, 1, radius, gram, histEdge, nullptr) == 0) {
    return fail("rdf edges");
  }
  const std::string outdir = "var/run/dseams-lang-out";
  if (seams_frame_rdf2d(built, 3.0, 0.5, nullptr, radius, gram) != 0 ||
      seams_frame_rdf2d(built, 3.0, 0.5, outdir.c_str(), radius, gram) != 0 ||
      seams_frame_rdf2d(nullptr, 3.0, 0.5, nullptr, radius, gram) == 0) {
    return fail("rdf2d");
  }
  int nrows = 0;
  if (seams_frame_prisms(built, 3.0, 1, 6, nullptr, nullptr, &nrows) != 0 || nrows < 1) {
    return fail("prism count");
  }
  std::vector<seams_prism_row> prows(static_cast<std::size_t>(nrows));
  if (seams_frame_prisms(nullptr, 3.0, 1, 2, nullptr, nullptr, nullptr) == 0 ||
      seams_frame_prisms(built, 3.0, 1, 6, outdir.c_str(), prows.data(), &nrows) != 0) {
    return fail("prisms");
  }
  int mono[9] = {};
  if (seams_frame_monolayer(built, 3.0, 1, 6, 10.0, nullptr, mono, 9) != 0 ||
      seams_frame_monolayer(built, 3.0, 1, 6, 10.0, outdir.c_str(), mono, 9) != 0 ||
      seams_frame_monolayer(nullptr, 3.0, 1, 6, 10.0, nullptr, mono, 9) == 0) {
    return fail("monolayer");
  }
  int nHex = 0, nCub = 0, nWat = 0, nBoth = 0;
  int perAtom[4] = {};
  if (seams_frame_cages(built, 3.0, 1, 4, "cutoff", 1, perAtom, &nHex, &nCub, &nWat, &nBoth) != 0 ||
      seams_frame_cages(nullptr, 3.0, 1, 4, "seeded", 0, nullptr, &nHex, &nCub, &nWat, &nBoth) == 0) {
    return fail("cage edges");
  }
  double qlEdge[4] = {};
  double qlbarEdge[4] = {};
  if (seams_frame_steinhardt(nullptr, 3.0, 1, 4, qlEdge, qlbarEdge) == 0 ||
      seams_frame_steinhardt(built, 3.0, 1, 4, nullptr, nullptr) != 0) {
    return fail("steinhardt edges");
  }
  int nCages = 0, nAtoms = 0;
  if (seams_frame_signature(nullptr, 3.0, 1, "cutoff", "512", &nCages, &nAtoms) == 0 ||
      seams_frame_signature(built, 3.0, 1, "nope", "512", &nCages, &nAtoms) == 0) {
    return fail("signature edges");
  }
  double tiltedBox[9] = {};
  seams_frame_box(tilted, tiltedBox);
  if (tiltedBox[6] == 0.0) {
    return fail("shear was dropped");
  }
  if (seams_path_nframes("nodot") < 0 && seams_last_error() == nullptr) {
    return fail("extension");
  }
  int guestHex = 0, guestCub = 0, guestWat = 0, guestBoth = 0;
  if (seams_frame_cages(guest, 3.5, 1, 0, "cutoff", 0, nullptr, &guestHex,
                        &guestCub, &guestWat, &guestBoth) != 0) {
    return fail("guest cutoff cages");
  }
  setenv("SEAMS_OFFLOAD", "1", 1);
  if (seams_frame_cages(guest, 3.5, 1, 4, "mutual", 0, nullptr, &guestHex,
                        &guestCub, &guestWat, &guestBoth) != 0) {
    return fail("guest offload cages");
  }
  unsetenv("SEAMS_OFFLOAD");
  const std::string emptyDump = "var/run/dseams-lang-empty.lammpstrj";
  {
    std::ofstream out(emptyDump);
    out << "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n0\n"
           "ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
           "ITEM: ATOMS id type x y z\n";
  }
  seams_frame *empty = seams_frame_read(emptyDump.c_str(), 1, -1);
  if (empty == nullptr || seams_frame_cages(empty, 3.0, 0, 4, "cutoff", 0, nullptr,
                                            &guestHex, &guestCub, &guestWat,
                                            &guestBoth) != 0) {
    return fail("empty cages");
  }
  const std::string tinyXyz = "var/run/dseams-lang-tiny.xyz";
  {
    std::ofstream out(tinyXyz);
    out << "1\ncomment\nO 0 0 0\n";
  }
  seams_frame *fromXyz = seams_frame_read(tinyXyz.c_str(), 1, -1);
  if (fromXyz == nullptr) {
    return fail("xyz read");
  }
  if (seams_frame_read("var/run/dseams-lang-missing.xyz", 1, -1) != nullptr) {
    return fail("missing xyz should fail");
  }
  if (seams_frame_chill(tip, 3.5, 99, "cutoff", labels.data(), counts) != 0) {
    return fail("chill type filter");
  }
  if (seams_frame_cages(built, 3.0, 99, 4, "nope", 0, nullptr, &guestHex, &guestCub,
                        &guestWat, &guestBoth) == 0) {
    return fail("cages should reject the graph");
  }
  if (seams_frame_cages(built, 3.0, 99, 4, "cutoff", 0, nullptr, &guestHex, &guestCub,
                        &guestWat, &guestBoth) != 0) {
    return fail("cages type filter");
  }
  seams_prism_row oneRow{};
  int shortRows = 1;
  if (seams_frame_prisms(built, 3.0, 1, 6, nullptr, &oneRow, &shortRows) == 0) {
    return fail("short prism buffer");
  }
  if (seams_frame_steinhardt(built, 3.0, 1, -1, qlEdge, qlbarEdge) != 0 &&
      seams_last_error() == nullptr) {
    return fail("steinhardt error");
  }
  if (seams_frame_rdf(built, -1.0, 4, 1, 1, radius, gram, histEdge, nullptr) != 0 &&
      seams_last_error() == nullptr) {
    return fail("rdf error");
  }
  if (seams_frame_rdf2d(built, -1.0, 0.1, nullptr, radius, gram) != 0 &&
      seams_last_error() == nullptr) {
    return fail("rdf2d error");
  }
  if (seams_frame_monolayer(built, -1.0, 1, 6, 1.0, nullptr, mono, 9) != 0 &&
      seams_last_error() == nullptr) {
    return fail("monolayer error");
  }
  seams_frame_free(fromXyz);
  seams_frame_free(empty);
  seams_frame_free(built);
  seams_frame_free(tilted);
  seams_frame_free(nullptr);

  seams_frame_free(molten);
  seams_frame_free(tip);
  seams_frame_free(ice);
  seams_frame_free(guest);
  std::cout << "seams_lang ok " << seams_lang_version() << "\n";
  return 0;
}
