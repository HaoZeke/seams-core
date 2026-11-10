// Copyright (c) d-SEAMS developers
#include <seams_insitu.h>

#include <bop.hpp>
#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <memory>
#ifdef SEAMS_HAS_OPENMP
#include <omp.h>
#endif
#include <sstream>
#include <string>
#include <vector>

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

const char *kChillName[SEAMS_INSITU_CHILL_N] = {
    "cubic",     "hexagonal", "water",       "interfacial", "clathrate",
    "interClathrate", "unclassified", "reCubic", "reHex"};

struct Engine {
  seams_insitu_options opt;
  seams_insitu_sink sink;
  void *sink_user;
  FILE *record;
  int last_status;
  int step;
  int n;
  int n_rings;
  int rings[SEAMS_INSITU_RING_MAX + 1];
  int chill[SEAMS_INSITU_CHILL_N];
  int n_hc;
  int n_ddc;
  int n_both;
  int n_water;
  int n_water_atoms;
  std::vector<int> labels;
  std::vector<int> cages;
  std::string text;
  seams_insitu_result result;
  bool have_result;
  // Reused when the next frame's bonds are almost the previous frame's.
  std::unique_ptr<primitive::RingUpdater> rings_up;
  std::unique_ptr<ring::AffiliationUpdater> cages_up;
  int rings_depth = -1;
  // Neighbours out to cutoff + skin. A later frame filters this list
  // while every atom has moved less than half the skin.
  std::vector<std::vector<int>> wide;
  std::vector<double> built_xyz;
  double built_lx = 0;
  double built_ly = 0;
  double built_lz = 0;
  std::vector<int> host_counts;
  std::vector<int> host_flat;
  bool host_fresh = false;
};

enum class Gate { Proceed, Done };

int fail(Engine *eng, int code) {
  eng->last_status = code;
  return code;
}

bool cadence_ok(const seams_insitu_cadence &c) {
  if (c.every < 1 || c.first < 0 || c.last < -1) {
    return false;
  }
  return true;
}

bool options_ok(const seams_insitu_options &opt) {
  if (!(opt.cutoff > 0.0)) {
    return false;
  }
  if (opt.ring_depth < 3 || opt.ring_depth > SEAMS_INSITU_RING_MAX) {
    return false;
  }
  if (opt.write_labels != 0 && opt.write_labels != 1) {
    return false;
  }
  if (opt.debug_dump != 0 && opt.debug_dump != 1) {
    return false;
  }
  return cadence_ok(opt.cadence);
}

std::string fmt17(double value) {
  char buf[64];
  const auto result =
      std::to_chars(buf, buf + sizeof(buf), value, std::chars_format::general, 17);
  return std::string(buf, result.ptr);
}

bool box_of(const seams_insitu_frame *frame, std::vector<double> &box,
            std::vector<double> &boxLow) {
  return nneigh::setBoxFromLammpsEdges(box, boxLow, frame->xlo, frame->ylo,
                                      frame->zlo, frame->lx, frame->ly,
                                      frame->lz, frame->xy, frame->xz,
                                      frame->yz);
}

int water_type_of(const Engine *eng, const Cloud &cloud) {
  if (eng->opt.type > 0) {
    return eng->opt.type;
  }
  if (cloud.nop > 0) {
    return cloud.pts[0].type;
  }
  return 1;
}

void publish(Engine *eng, int status) {
  eng->last_status = status;
  eng->have_result = true;
  eng->result.status = status;
  eng->result.step = eng->step;
  eng->result.n = eng->n;
  eng->result.n_rings = eng->n_rings;
  for (int i = 0; i <= SEAMS_INSITU_RING_MAX; ++i) {
    eng->result.ring_count[i] = eng->rings[i];
  }
  for (int i = 0; i < SEAMS_INSITU_CHILL_N; ++i) {
    eng->result.chill_count[i] = eng->chill[i];
  }
  eng->result.n_hc = eng->n_hc;
  eng->result.n_ddc = eng->n_ddc;
  eng->result.n_both = eng->n_both;
  eng->result.n_water = eng->n_water;
  eng->result.n_water_atoms = eng->n_water_atoms;
  eng->result.chill_labels = eng->labels.empty() ? nullptr : eng->labels.data();
  eng->result.cage_codes = eng->cages.empty() ? nullptr : eng->cages.data();
  eng->result.text = eng->text.c_str();
  eng->result.text_len = static_cast<int>(eng->text.size());
}

int notify(Engine *eng, int status) {
  publish(eng, status);
  if (eng->sink == nullptr) {
    return status;
  }
  const int rc = eng->sink(&eng->result, eng->sink_user);
  if (rc != 0) {
    eng->last_status = SEAMS_INSITU_ERR_IO;
    eng->result.status = SEAMS_INSITU_ERR_IO;
    return SEAMS_INSITU_ERR_IO;
  }
  return status;
}

void compose(Engine *eng) {
  std::ostringstream os;
  os << "step " << eng->step << "\n";
  os << "n " << eng->n << "\n";
  os << "rings";
  for (int size = 3; size <= eng->opt.ring_depth; ++size) {
    os << " " << size << ":" << eng->rings[size];
  }
  os << "\n";
  os << "n_rings " << eng->n_rings << "\n";
  os << "chill";
  for (int i = 0; i < SEAMS_INSITU_CHILL_N; ++i) {
    os << " " << kChillName[i] << ":" << eng->chill[i];
  }
  os << "\n";
  os << "cage hc:" << eng->n_hc << " ddc:" << eng->n_ddc
     << " both:" << eng->n_both << " water:" << eng->n_water
     << " n_water:" << eng->n_water_atoms << "\n";
  if (eng->opt.write_labels) {
    os << "labels";
    for (int label : eng->labels) {
      os << " " << label;
    }
    os << "\n";
    os << "cages";
    for (int code : eng->cages) {
      os << " " << code;
    }
    os << "\n";
  }
  eng->text = os.str();
}

void clear_counts(Engine *eng) {
  eng->n_rings = 0;
  for (int i = 0; i <= SEAMS_INSITU_RING_MAX; ++i) {
    eng->rings[i] = 0;
  }
  for (int i = 0; i < SEAMS_INSITU_CHILL_N; ++i) {
    eng->chill[i] = 0;
  }
  eng->n_hc = eng->n_ddc = eng->n_both = eng->n_water = eng->n_water_atoms = 0;
  eng->labels.clear();
  eng->cages.clear();
}

// Pairs in a wide id-list whose current minimum image is inside the cutoff.
// The comparison matches neighListO: a distance equal to the cutoff is out.
std::vector<std::vector<int>> filter_cutoff(const Cloud &cloud,
                                            const std::vector<std::vector<int>> &wide,
                                            double cutoff, int water) {
  const double limit = cutoff * cutoff;
  const gen::FracBox box = gen::makeFracBox(cloud);
  int maxId = 0;
  for (int i = 0; i < cloud.nop; ++i) {
    maxId = std::max(maxId, cloud.pts[static_cast<std::size_t>(i)].atomID);
  }
  std::vector<int> at(static_cast<std::size_t>(maxId) + 1, -1);
  for (int i = 0; i < cloud.nop; ++i) {
    const int id = cloud.pts[static_cast<std::size_t>(i)].atomID;
    if (id >= 0) {
      at[static_cast<std::size_t>(id)] = i;
    }
  }
  const auto indexOf = [&](int id) {
    return (id >= 0 && id <= maxId) ? at[static_cast<std::size_t>(id)] : -1;
  };
  static const bool phase = std::getenv("DSEAMS_PHASE") != nullptr;
  const auto filtT = phase ? std::chrono::steady_clock::now()
                           : std::chrono::steady_clock::time_point{};
  std::vector<std::vector<int>> out(static_cast<std::size_t>(cloud.nop));
  for (int i = 0; i < cloud.nop; ++i) {
    out[static_cast<std::size_t>(i)].push_back(cloud.pts[static_cast<std::size_t>(i)].atomID);
  }
  const int nRows = static_cast<int>(wide.size());
  for (int r = 0; r < nRows; ++r) {
    if (!wide[static_cast<std::size_t>(r)].empty()) {
      const int iatom = indexOf(wide[static_cast<std::size_t>(r)][0]);
      if (iatom >= 0) {
        out[static_cast<std::size_t>(iatom)].reserve(
            wide[static_cast<std::size_t>(r)].size());
      }
    }
  }
#ifdef SEAMS_HAS_OPENMP
#pragma omp parallel for schedule(static) if (nRows >= 256 && !omp_in_parallel())
#endif
  for (int r = 0; r < nRows; ++r) {
    const auto &row = wide[static_cast<std::size_t>(r)];
    if (row.empty()) {
      continue;
    }
    const int iatom = indexOf(row[0]);
    if (iatom < 0 || cloud.pts[static_cast<std::size_t>(iatom)].type != water) {
      continue;
    }
    auto &dest = out[static_cast<std::size_t>(iatom)];
    for (std::size_t k = 1; k < row.size(); ++k) {
      const int jatom = indexOf(row[k]);
      if (jatom < 0 || jatom == iatom ||
          cloud.pts[static_cast<std::size_t>(jatom)].type != water) {
        continue;
      }
      if (gen::periodicDistSq(box, cloud, iatom, jatom) >= limit) {
        continue;
      }
      // A host list can name one atom through several periodic images.
      if (std::find(dest.begin() + 1, dest.end(), row[k]) != dest.end()) {
        continue;
      }
      dest.push_back(row[k]);
    }
  }
  if (phase) {
    std::fprintf(stderr, "phase filt %.6f\n",
                 std::chrono::duration<double, std::milli>(
                     std::chrono::steady_clock::now() - filtT)
                     .count());
  }
  return out;
}

std::vector<std::vector<int>> rows_from_flat(const Cloud &cloud, const std::vector<int> &counts,
                                            const std::vector<int> &flat) {
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(cloud.nop));
  std::size_t cursor = 0;
  for (const int count : counts) {
    const int self = flat[cursor];
    const auto at = cloud.idIndexMap.find(self);
    if (at != cloud.idIndexMap.end()) {
      auto &row = rows[static_cast<std::size_t>(at->second)];
      row.clear();
      const std::size_t n = std::min(static_cast<std::size_t>(count), flat.size() - cursor);
      row.insert(row.end(), flat.begin() + static_cast<std::ptrdiff_t>(cursor),
                 flat.begin() + static_cast<std::ptrdiff_t>(cursor + n));
    }
    cursor += static_cast<std::size_t>(count);
  }
  return rows;
}

void remember_build(Engine *eng, const Cloud &cloud, std::vector<std::vector<int>> wide) {
  eng->wide = std::move(wide);
  eng->built_xyz.resize(static_cast<std::size_t>(cloud.nop) * 3);
  for (int i = 0; i < cloud.nop; ++i) {
    const auto &p = cloud.pts[static_cast<std::size_t>(i)];
    eng->built_xyz[static_cast<std::size_t>(3 * i)] = p.x;
    eng->built_xyz[static_cast<std::size_t>(3 * i + 1)] = p.y;
    eng->built_xyz[static_cast<std::size_t>(3 * i + 2)] = p.z;
  }
  eng->built_lx = cloud.box[0];
  eng->built_ly = cloud.box[1];
  eng->built_lz = cloud.box[2];
}

int analyze_cloud(Engine *eng, Cloud &cloud, int step, bool do_sink) {
  eng->step = step;
  eng->n = cloud.nop;
  clear_counts(eng);
  const int water = water_type_of(eng, cloud);
  if (cloud.nop > 0) {
    try {
      static const bool phase = std::getenv("DSEAMS_PHASE") != nullptr;
      const auto phaseStart = std::chrono::steady_clock::now();
      auto stamp = phaseStart;
      auto mark = [&](const char *name) {
        if (!phase) {
          return;
        }
        const auto now = std::chrono::steady_clock::now();
        std::fprintf(stderr, "phase %s %.6f\n", name,
                     std::chrono::duration<double, std::milli>(now - stamp).count());
        stamp = now;
      };
      std::vector<std::vector<int>> nList;
      const double skin = 0.5;
      const double half = skin * 0.5;
      bool keep = cloud.box.size() >= 3 &&
                  eng->wide.size() == static_cast<std::size_t>(cloud.nop) &&
                  eng->built_lx == cloud.box[0] && eng->built_ly == cloud.box[1] &&
                  eng->built_lz == cloud.box[2];
      if (keep) {
        for (int i = 0; i < cloud.nop; ++i) {
          const auto &p = cloud.pts[static_cast<std::size_t>(i)];
          const double dx = p.x - eng->built_xyz[static_cast<std::size_t>(3 * i)];
          const double dy = p.y - eng->built_xyz[static_cast<std::size_t>(3 * i + 1)];
          const double dz = p.z - eng->built_xyz[static_cast<std::size_t>(3 * i + 2)];
          if (dx * dx + dy * dy + dz * dz > half * half) {
            keep = false;
            break;
          }
        }
      }
      if (keep) {
        nList = filter_cutoff(cloud, eng->wide, eng->opt.cutoff, water);
        eng->host_fresh = false;
      } else if (eng->host_fresh) {
        const auto shell = filter_cutoff(
            cloud, rows_from_flat(cloud, eng->host_counts, eng->host_flat),
            eng->opt.cutoff + skin, water);
        nList = filter_cutoff(cloud, shell, eng->opt.cutoff, water);
        eng->host_fresh = false;
        remember_build(eng, cloud, shell);
      } else {
        auto wide = nneigh::neighListO(eng->opt.cutoff + skin, cloud, water);
        nList = filter_cutoff(cloud, wide, eng->opt.cutoff, water);
        remember_build(eng, cloud, std::move(wide));
      }
      mark("neigh");
      auto idx = nneigh::neighbourListByIndex(cloud, nList);
      if (eng->rings_up == nullptr) {
        eng->rings_up = std::make_unique<primitive::RingUpdater>(eng->opt.ring_depth);
        eng->rings_up->setFullPassFraction(0.15);
        eng->rings_depth = eng->opt.ring_depth;
      }
      const auto &rings = eng->rings_up->update(idx);
      mark("rings");
      for (const auto &ring : rings) {
        const int size = static_cast<int>(ring.size());
        eng->rings[size] += 1;
        eng->n_rings += 1;
      }
      std::vector<std::vector<int>> six;
      for (const auto &ring : rings) {
        if (ring.size() == 6) {
          six.push_back(ring);
        }
      }
      mark("six");
      const auto faces = ring::contractibleRings(six, cloud);
      mark("contract");
      if (eng->cages_up == nullptr) {
        eng->cages_up = std::make_unique<ring::AffiliationUpdater>();
        eng->cages_up->setBatchFraction(0.5);
      }
      const auto &aff = eng->cages_up->update(faces, idx);
      mark("affil");
      std::vector<int> hc(static_cast<std::size_t>(cloud.nop), 0);
      std::vector<int> ddc(static_cast<std::size_t>(cloud.nop), 0);
      for (std::size_t r = 0; r < faces.size(); ++r) {
        const int h = (r < aff.hc.size() && aff.hc[r]) ? 1 : 0;
        const int d = (r < aff.ddc.size() && aff.ddc[r]) ? 1 : 0;
        for (const int atom : faces[r]) {
          if (atom >= 0 && atom < cloud.nop) {
            const std::size_t a = static_cast<std::size_t>(atom);
            hc[a] = hc[a] | h;
            ddc[a] = ddc[a] | d;
          }
        }
      }
      mark("scatter");
      chill::getCorrelPlus(cloud, nList, false);
      mark("correl");
      chill::getIceTypePlusNoPrint(cloud, nList, false);
      mark("icetype");
      chill::clearIceTypesExcept(cloud, std::vector<int>{water});
      mark("clear");
      eng->labels.resize(static_cast<std::size_t>(cloud.nop));
      eng->cages.assign(static_cast<std::size_t>(cloud.nop), 0);
      std::vector<int> types(static_cast<std::size_t>(cloud.nop));
      for (int i = 0; i < cloud.nop; ++i) {
        const int label = static_cast<int>(cloud.pts[static_cast<std::size_t>(i)].iceType);
        eng->labels[static_cast<std::size_t>(i)] = label;
        eng->chill[label] += 1;
        types[static_cast<std::size_t>(i)] = cloud.pts[static_cast<std::size_t>(i)].type;
      }
      seams_insitu_tally_cages(cloud.nop, types.data(), water, hc.data(), ddc.data(),
                               eng->cages.data(), &eng->n_hc, &eng->n_ddc,
                               &eng->n_both, &eng->n_water);
      eng->n_water_atoms = 0;
      for (int i = 0; i < cloud.nop; ++i) {
        if (cloud.pts[static_cast<std::size_t>(i)].type == water) {
          eng->n_water_atoms += 1;
        }
      }
    } catch (const std::exception &) {
      return fail(eng, SEAMS_INSITU_ERR);
    }
  }
  compose(eng);
  if (!do_sink) {
    publish(eng, SEAMS_INSITU_OK);
    return SEAMS_INSITU_OK;
  }
  return notify(eng, SEAMS_INSITU_OK);
}

int build_cloud(const Engine *eng, const seams_insitu_frame *frame, Cloud &cloud) {
  cloud = Cloud{};
  cloud.nop = frame->n;
  cloud.currentFrame = frame->step;
  cloud.pts.resize(static_cast<std::size_t>(frame->n));
  const double scale = frame->coord_scale > 0.0 ? frame->coord_scale : 1.0;
  for (int i = 0; i < frame->n; ++i) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.x = frame->xyz[3 * i] * scale;
    point.y = frame->xyz[3 * i + 1] * scale;
    point.z = frame->xyz[3 * i + 2] * scale;
    point.type = frame->types != nullptr ? frame->types[i] : 1;
    point.atomID = frame->ids != nullptr ? frame->ids[i] : i + 1;
    point.molID = point.atomID;
    point.inSlice = true;
    if (cloud.idIndexMap.find(point.atomID) != cloud.idIndexMap.end()) {
      return SEAMS_INSITU_ERR;
    }
    cloud.idIndexMap[point.atomID] = i;
  }
  box_of(frame, cloud.box, cloud.boxLow);
  (void)eng;
  return SEAMS_INSITU_OK;
}

Gate gate(Engine *eng, const seams_insitu_frame *frame, int *status) {
  if (eng == nullptr || frame == nullptr) {
    if (eng != nullptr) {
      *status = fail(eng, SEAMS_INSITU_ERR);
    } else {
      *status = SEAMS_INSITU_ERR;
    }
    return Gate::Done;
  }
  if (!seams_insitu_due(&eng->opt.cadence, frame->step)) {
    *status = fail(eng, SEAMS_INSITU_SKIP);
    return Gate::Done;
  }
  if (frame->gpu_resident) {
    eng->step = frame->step;
    eng->n = frame->n;
    clear_counts(eng);
    std::ostringstream os;
    os << "step " << frame->step << "\n";
    os << "gpu deferred\n";
    os << "device " << frame->gpu_device << "\n";
    os << "n " << frame->n << "\n";
    eng->text = os.str();
    *status = notify(eng, SEAMS_INSITU_GPU_DEFERRED);
    return Gate::Done;
  }
  if (frame->local_only && !frame->gathered) {
    *status = fail(eng, SEAMS_INSITU_NEED_GATHER);
    return Gate::Done;
  }
  if (frame->n < 0 || (frame->n > 0 && frame->xyz == nullptr)) {
    *status = fail(eng, SEAMS_INSITU_ERR);
    return Gate::Done;
  }
  if (frame->n > 0) {
    std::vector<double> box;
    std::vector<double> boxLow;
    if (!box_of(frame, box, boxLow)) {
      *status = fail(eng, SEAMS_INSITU_ERR_BOX);
      return Gate::Done;
    }
  }
  return Gate::Proceed;
}

int analyze_frame(Engine *eng, const seams_insitu_frame *frame) {
  if (frame->n == 0) {
    Cloud empty;
    return analyze_cloud(eng, empty, frame->step, true);
  }
  Cloud cloud;
  const int built = build_cloud(eng, frame, cloud);
  if (built != SEAMS_INSITU_OK) {
    return fail(eng, built);
  }
  return analyze_cloud(eng, cloud, frame->step, true);
}

bool write_frame(std::ostream &out, const seams_insitu_frame *frame,
                 const std::vector<double> &box, const std::vector<double> &boxLow) {
  const bool tilt = box.size() >= 6;
  out << "ITEM: TIMESTEP\n" << frame->step << "\n";
  out << "ITEM: NUMBER OF ATOMS\n" << frame->n << "\n";
  if (tilt) {
    out << "ITEM: BOX BOUNDS xy xz yz pp pp pp\n";
  } else {
    out << "ITEM: BOX BOUNDS pp pp pp\n";
  }
  for (int axis = 0; axis < 3; ++axis) {
    const double lo = boxLow[static_cast<std::size_t>(axis)];
    const double hi = lo + box[static_cast<std::size_t>(axis)];
    out << fmt17(lo) << " " << fmt17(hi);
    if (tilt) {
      out << " " << fmt17(box[static_cast<std::size_t>(3 + axis)]);
    }
    out << "\n";
  }
  out << "ITEM: ATOMS id type x y z\n";
  for (int i = 0; i < frame->n; ++i) {
    const int id = frame->ids != nullptr ? frame->ids[i] : i + 1;
    const int type = frame->types != nullptr ? frame->types[i] : 1;
    out << id << " " << type << " " << fmt17(frame->xyz[3 * i]) << " "
        << fmt17(frame->xyz[3 * i + 1]) << " " << fmt17(frame->xyz[3 * i + 2])
        << "\n";
  }
  return static_cast<bool>(out);
}

std::vector<int> timesteps_of(const std::string &path, bool *ok) {
  std::ifstream in(path);
  std::vector<int> steps;
  if (!in) {
    *ok = false;
    return steps;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.compare(0, 14, "ITEM: TIMESTEP") != 0) {
      continue;
    }
    if (!std::getline(in, line)) {
      *ok = false;
      return steps;
    }
    int step = 0;
    const auto parsed = std::from_chars(line.data(), line.data() + line.size(), step);
    if (parsed.ec != std::errc()) {
      *ok = false;
      return steps;
    }
    steps.push_back(step);
  }
  *ok = true;
  return steps;
}

int file_sink(const seams_insitu_result *result, void *user) {
  auto *file = static_cast<FILE *>(user);
  if (file == nullptr || result == nullptr || result->text == nullptr) {
    return 1;
  }
  const auto wrote = std::fwrite(result->text, 1, static_cast<std::size_t>(result->text_len), file);
  if (wrote != static_cast<std::size_t>(result->text_len)) {
    return 1;
  }
  if (std::fflush(file) != 0) {
    return 1;
  }
  return 0;
}

} // namespace

void seams_insitu_options_default(seams_insitu_options *opt) {
  if (opt == nullptr) {
    return;
  }
  opt->cutoff = 3.5;
  opt->type = 1;
  opt->ring_depth = 8;
  opt->write_labels = 1;
  opt->debug_dump = 0;
  opt->cadence.every = 1;
  opt->cadence.first = 0;
  opt->cadence.last = -1;
}

seams_insitu_engine *seams_insitu_create(const seams_insitu_options *opt) {
  seams_insitu_options local;
  seams_insitu_options_default(&local);
  if (opt != nullptr) {
    local = *opt;
  }
  if (!options_ok(local)) {
    return nullptr;
  }
  auto *eng = new Engine();
  eng->opt = local;
  eng->sink = nullptr;
  eng->sink_user = nullptr;
  eng->record = nullptr;
  eng->last_status = SEAMS_INSITU_OK;
  eng->step = 0;
  eng->n = 0;
  eng->have_result = false;
  clear_counts(eng);
  std::memset(&eng->result, 0, sizeof(eng->result));
  return reinterpret_cast<seams_insitu_engine *>(eng);
}

void seams_insitu_destroy(seams_insitu_engine *eng) {
  if (eng == nullptr) {
    return;
  }
  auto *impl = reinterpret_cast<Engine *>(eng);
  if (impl->record != nullptr) {
    std::fclose(impl->record);
    impl->record = nullptr;
  }
  delete impl;
}

int seams_insitu_set_sink(seams_insitu_engine *eng, seams_insitu_sink sink, void *user) {
  if (eng == nullptr) {
    return SEAMS_INSITU_ERR;
  }
  auto *impl = reinterpret_cast<Engine *>(eng);
  impl->sink = sink;
  impl->sink_user = user;
  return SEAMS_INSITU_OK;
}

int seams_insitu_open_record(seams_insitu_engine *eng, const char *path) {
  if (eng == nullptr || path == nullptr || path[0] == '\0') {
    return SEAMS_INSITU_ERR;
  }
  auto *impl = reinterpret_cast<Engine *>(eng);
  FILE *file = std::fopen(path, "w");
  if (file == nullptr) {
    return fail(impl, SEAMS_INSITU_ERR_IO);
  }
  if (impl->record != nullptr) {
    std::fclose(impl->record);
  }
  impl->record = file;
  impl->sink = file_sink;
  impl->sink_user = file;
  return SEAMS_INSITU_OK;
}

int seams_insitu_close_record(seams_insitu_engine *eng) {
  if (eng == nullptr) {
    return SEAMS_INSITU_ERR;
  }
  auto *impl = reinterpret_cast<Engine *>(eng);
  if (impl->record == nullptr) {
    return SEAMS_INSITU_OK;
  }
  std::fclose(impl->record);
  impl->record = nullptr;
  impl->sink = nullptr;
  impl->sink_user = nullptr;
  return SEAMS_INSITU_OK;
}

int seams_insitu_due(const seams_insitu_cadence *cadence, int step) {
  if (cadence == nullptr || !cadence_ok(*cadence)) {
    return -1;
  }
  if (step < cadence->first) {
    return 0;
  }
  if (cadence->last >= 0 && step > cadence->last) {
    return 0;
  }
  if (((step - cadence->first) % cadence->every) != 0) {
    return 0;
  }
  return 1;
}

int seams_insitu_append_dump(const seams_insitu_frame *frame, const char *path) {
  if (frame == nullptr || path == nullptr || path[0] == '\0') {
    return SEAMS_INSITU_ERR;
  }
  if (frame->n < 0 || (frame->n > 0 && frame->xyz == nullptr)) {
    return SEAMS_INSITU_ERR;
  }
  std::vector<double> box;
  std::vector<double> boxLow;
  if (frame->n > 0 && !box_of(frame, box, boxLow)) {
    return SEAMS_INSITU_ERR_BOX;
  }
  if (frame->n == 0) {
    box = {1.0, 1.0, 1.0};
    boxLow = {0.0, 0.0, 0.0};
  }
  std::ofstream out(path, std::ios::app);
  if (!out) {
    return SEAMS_INSITU_ERR_IO;
  }
  if (!write_frame(out, frame, box, boxLow)) {
    return SEAMS_INSITU_ERR_IO;
  }
  return SEAMS_INSITU_OK;
}

int seams_insitu_set_neighbours(seams_insitu_engine *eng, int n, const int *counts,
                                const int *flat) {
  auto *impl = reinterpret_cast<Engine *>(eng);
  if (impl == nullptr || n < 0) {
    return SEAMS_INSITU_ERR;
  }
  if (n == 0) {
    impl->host_counts.clear();
    impl->host_flat.clear();
    impl->host_fresh = false;
    impl->wide.clear();
    impl->built_xyz.clear();
    return SEAMS_INSITU_OK;
  }
  if (counts == nullptr || flat == nullptr) {
    return SEAMS_INSITU_ERR;
  }
  std::size_t total = 0;
  for (int i = 0; i < n; ++i) {
    if (counts[i] < 1) {
      return SEAMS_INSITU_ERR;
    }
    total += static_cast<std::size_t>(counts[i]);
  }
  impl->host_counts.assign(counts, counts + n);
  impl->host_flat.assign(flat, flat + total);
  impl->host_fresh = true;
  return SEAMS_INSITU_OK;
}

int seams_insitu_step(seams_insitu_engine *eng, const seams_insitu_frame *frame) {
  auto *impl = reinterpret_cast<Engine *>(eng);
  int status = SEAMS_INSITU_OK;
  if (gate(impl, frame, &status) == Gate::Done) {
    return status;
  }
  return analyze_frame(impl, frame);
}

int seams_insitu_observe(seams_insitu_engine *eng, const seams_insitu_frame *frame,
                         const char *dump_path) {
  auto *impl = reinterpret_cast<Engine *>(eng);
  int status = SEAMS_INSITU_OK;
  if (gate(impl, frame, &status) == Gate::Done) {
    return status;
  }
  const int rc = analyze_frame(impl, frame);
  if (rc != SEAMS_INSITU_OK) {
    return rc;
  }
  if (impl->opt.debug_dump == 0 || dump_path == nullptr || dump_path[0] == '\0') {
    return SEAMS_INSITU_OK;
  }
  const int wrote = seams_insitu_append_dump(frame, dump_path);
  if (wrote != SEAMS_INSITU_OK) {
    return fail(impl, wrote);
  }
  return SEAMS_INSITU_OK;
}

int seams_insitu_replay(seams_insitu_engine *eng, const char *dump_path,
                        const char *record_path) {
  if (eng == nullptr || dump_path == nullptr || record_path == nullptr) {
    return SEAMS_INSITU_ERR;
  }
  auto *impl = reinterpret_cast<Engine *>(eng);
  bool steps_ok = false;
  const auto steps = timesteps_of(dump_path, &steps_ok);
  if (!steps_ok) {
    return fail(impl, SEAMS_INSITU_ERR_IO);
  }
  if (steps.empty()) {
    return fail(impl, SEAMS_INSITU_ERR);
  }
  sinp::dropLammpsDumpIndex(dump_path);
  std::ofstream out(record_path, std::ios::trunc);
  if (!out) {
    return fail(impl, SEAMS_INSITU_ERR_IO);
  }
  for (int frame = 1; frame <= static_cast<int>(steps.size()); ++frame) {
    const int step = steps[static_cast<std::size_t>(frame - 1)];
    if (!seams_insitu_due(&impl->opt.cadence, step)) {
      continue;
    }
    Cloud cloud;
    sinp::readLammpsTrj(dump_path, frame, cloud, false, {}, {});
    const int rc = analyze_cloud(impl, cloud, step, false);
    if (rc != SEAMS_INSITU_OK) {
      return rc;
    }
    out << impl->text;
    if (!out) {
      return fail(impl, SEAMS_INSITU_ERR_IO);
    }
  }
  return SEAMS_INSITU_OK;
}

int seams_insitu_status(const seams_insitu_engine *eng) {
  if (eng == nullptr) {
    return SEAMS_INSITU_ERR;
  }
  return reinterpret_cast<const Engine *>(eng)->last_status;
}

const seams_insitu_result *seams_insitu_last(const seams_insitu_engine *eng) {
  if (eng == nullptr) {
    return nullptr;
  }
  const auto *impl = reinterpret_cast<const Engine *>(eng);
  if (!impl->have_result) {
    return nullptr;
  }
  return &impl->result;
}

int seams_insitu_frame_bytes(void) { return static_cast<int>(sizeof(seams_insitu_frame)); }

int seams_insitu_options_bytes(void) {
  return static_cast<int>(sizeof(seams_insitu_options));
}

int seams_insitu_result_bytes(void) {
  return static_cast<int>(sizeof(seams_insitu_result));
}

int seams_insitu_pack_owned(int nlocal, const int *mask, int groupbit,
                            const double *xyz, const int *types, const int *ids,
                            double *xyz_out, int *types_out, int *ids_out) {
  if (nlocal < 0) {
    return -1;
  }
  if (nlocal == 0) {
    return 0;
  }
  if (xyz == nullptr || xyz_out == nullptr || types_out == nullptr || ids_out == nullptr) {
    return -1;
  }
  int kept = 0;
  for (int i = 0; i < nlocal; ++i) {
    if (mask != nullptr && (mask[i] & groupbit) == 0) {
      continue;
    }
    xyz_out[3 * kept] = xyz[3 * i];
    xyz_out[3 * kept + 1] = xyz[3 * i + 1];
    xyz_out[3 * kept + 2] = xyz[3 * i + 2];
    types_out[kept] = types != nullptr ? types[i] : 1;
    ids_out[kept] = ids != nullptr ? ids[i] : i + 1;
    kept += 1;
  }
  return kept;
}

int seams_insitu_gather_displs(int nranks, const int *counts, int *displs, int *total) {
  if (nranks < 0 || total == nullptr) {
    return 1;
  }
  if (nranks > 0 && (counts == nullptr || displs == nullptr)) {
    return 1;
  }
  int sum = 0;
  for (int rank = 0; rank < nranks; ++rank) {
    if (counts[rank] < 0) {
      return 1;
    }
    displs[rank] = sum;
    sum += counts[rank];
  }
  *total = sum;
  return 0;
}

void seams_insitu_scale_nm(const double *xyz_nm, int n, double *xyz_a) {
  if (xyz_nm == nullptr || xyz_a == nullptr || n <= 0) {
    return;
  }
  for (int i = 0; i < 3 * n; ++i) {
    xyz_a[i] = xyz_nm[i] * 10.0;
  }
}

int seams_insitu_box_from_rows(const double rows[9], double *xlo, double *ylo, double *zlo,
                               double *lx, double *ly, double *lz, double *xy, double *xz,
                               double *yz) {
  if (rows == nullptr || xlo == nullptr || ylo == nullptr || zlo == nullptr || lx == nullptr ||
      ly == nullptr || lz == nullptr || xy == nullptr || xz == nullptr || yz == nullptr) {
    return 1;
  }
  const double ax = rows[0];
  const double ay = rows[1];
  const double az = rows[2];
  const double bx = rows[3];
  const double by = rows[4];
  const double bz = rows[5];
  const double cx = rows[6];
  const double cy = rows[7];
  const double cz = rows[8];
  const double scale = std::max(std::abs(ax), std::max(std::abs(by), std::abs(cz)));
  const double tol = 1e-8 * (scale > 1.0 ? scale : 1.0);
  if (!(ax > 0.0 && by > 0.0 && cz > 0.0)) {
    return 1;
  }
  if (std::abs(ay) > tol || std::abs(az) > tol || std::abs(bz) > tol) {
    return 1;
  }
  *xlo = 0.0;
  *ylo = 0.0;
  *zlo = 0.0;
  *lx = ax;
  *ly = by;
  *lz = cz;
  *xy = bx;
  *xz = cx;
  *yz = cy;
  return 0;
}

void seams_insitu_tally_cages(int n, const int *types, int water_type, const int *hc,
                              const int *ddc, int *codes, int *n_hc, int *n_ddc, int *n_both,
                              int *n_water) {
  if (n_hc == nullptr || n_ddc == nullptr || n_both == nullptr || n_water == nullptr) {
    return;
  }
  *n_hc = *n_ddc = *n_both = *n_water = 0;
  if (n <= 0 || types == nullptr || hc == nullptr || ddc == nullptr || codes == nullptr) {
    return;
  }
  for (int i = 0; i < n; ++i) {
    if (types[i] != water_type) {
      codes[i] = 0;
      continue;
    }
    const bool hexagonal = hc[i] != 0;
    const bool diamond = ddc[i] != 0;
    if (hexagonal && diamond) {
      codes[i] = 3;
      *n_both += 1;
    } else if (hexagonal) {
      codes[i] = 1;
      *n_hc += 1;
    } else if (diamond) {
      codes[i] = 2;
      *n_ddc += 1;
    } else {
      codes[i] = 0;
      *n_water += 1;
    }
  }
}

int seams_insitu_main(int argc, char **argv) {
  if (argc != 4 && argc != 5) {
    std::fprintf(stderr, "usage: seams-insitu replay DUMP RECORD [CUTOFF]\n");
    return 2;
  }
  if (std::strcmp(argv[1], "replay") != 0) {
    std::fprintf(stderr, "usage: seams-insitu replay DUMP RECORD [CUTOFF]\n");
    return 2;
  }
  seams_insitu_options opt;
  seams_insitu_options_default(&opt);
  if (argc == 5) {
    char *end = nullptr;
    opt.cutoff = std::strtod(argv[4], &end);
    if (end == argv[4] || *end != '\0') {
      std::fprintf(stderr, "cutoff is not a number\n");
      return 2;
    }
  }
  seams_insitu_engine *eng = seams_insitu_create(&opt);
  if (eng == nullptr) {
    std::fprintf(stderr, "rejected options\n");
    return 2;
  }
  const int rc = seams_insitu_replay(eng, argv[2], argv[3]);
  seams_insitu_destroy(eng);
  return rc == SEAMS_INSITU_OK ? 0 : 1;
}
