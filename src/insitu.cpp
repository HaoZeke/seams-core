#include <seams_insitu.h>

#include <bop.hpp>
#include <cage_affiliation.hpp>
#include <franzblau.hpp>
#include <mol_sys.hpp>
#include <neighbours.hpp>
#include <seams_input.hpp>
#include <topo_bulk.hpp>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
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

int analyze_cloud(Engine *eng, Cloud &cloud, int step, bool do_sink) {
  eng->step = step;
  eng->n = cloud.nop;
  clear_counts(eng);
  const int water = water_type_of(eng, cloud);
  if (cloud.nop > 0) {
    try {
      const auto nList = nneigh::neighListO(eng->opt.cutoff, cloud, water);
      const auto idx = nneigh::neighbourListByIndex(cloud, nList);
      const auto rings = primitive::ringNetwork(idx, eng->opt.ring_depth);
      for (const auto &ring : rings) {
        const int size = static_cast<int>(ring.size());
        if (size >= 0 && size <= SEAMS_INSITU_RING_MAX) {
          eng->rings[size] += 1;
        }
        eng->n_rings += 1;
      }
      std::vector<std::vector<int>> six;
      for (const auto &ring : rings) {
        if (ring.size() == 6) {
          six.push_back(ring);
        }
      }
      const auto faces = ring::contractibleRings(six, cloud);
      const auto aff = ring::cageAffiliation(faces, idx);
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
      chill::getCorrelPlus(cloud, nList, false);
      chill::getIceTypePlusNoPrint(cloud, nList, false);
      chill::clearIceTypesExcept(cloud, std::vector<int>{water});
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
  for (int i = 0; i < frame->n; ++i) {
    auto &point = cloud.pts[static_cast<std::size_t>(i)];
    point.x = frame->xyz[3 * i];
    point.y = frame->xyz[3 * i + 1];
    point.z = frame->xyz[3 * i + 2];
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
  if (dump_path == nullptr || dump_path[0] == '\0') {
    return analyze_frame(impl, frame);
  }
  const int wrote = seams_insitu_append_dump(frame, dump_path);
  if (wrote != SEAMS_INSITU_OK) {
    return fail(impl, wrote);
  }
  sinp::dropLammpsDumpIndex(dump_path);
  const int nframes = sinp::nLammpsFrames(dump_path);
  if (nframes < 1) {
    return fail(impl, SEAMS_INSITU_ERR_IO);
  }
  bool steps_ok = false;
  const auto steps = timesteps_of(dump_path, &steps_ok);
  if (!steps_ok || static_cast<int>(steps.size()) != nframes) {
    return fail(impl, SEAMS_INSITU_ERR);
  }
  Cloud cloud;
  sinp::readLammpsTrj(dump_path, nframes, cloud, false, {}, {});
  if (cloud.nop != frame->n) {
    return fail(impl, SEAMS_INSITU_ERR);
  }
  return analyze_cloud(impl, cloud, steps.back(), true);
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
