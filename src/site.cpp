//-----------------------------------------------------------------------------------
// d-SEAMS - Deferred Structural Elucidation Analysis for Molecular Simulations
// SPDX-License-Identifier: MIT
//-----------------------------------------------------------------------------------

#include <site.hpp>

#include <franzblau.hpp>
#include <generic.hpp>
#include <mol_sys.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace site {

namespace {

using Cloud = molSys::PointCloud<molSys::Point<double>, double>;

std::string trim(std::string_view s) {
  const auto begin = s.find_first_not_of(" \t\n\r");
  if (begin == std::string_view::npos) {
    return {};
  }
  const auto end = s.find_last_not_of(" \t\n\r");
  return std::string(s.substr(begin, end - begin + 1));
}

Kind kindFromName(std::string_view name) {
  if (name == "unspecified") {
    return Kind::unspecified;
  }
  if (name == "cationHead") {
    return Kind::cationHead;
  }
  if (name == "anion") {
    return Kind::anion;
  }
  if (name == "tail") {
    return Kind::tail;
  }
  if (name == "donorH") {
    return Kind::donorH;
  }
  if (name == "acceptor") {
    return Kind::acceptor;
  }
  if (name == "polar") {
    return Kind::polar;
  }
  if (name == "apolar") {
    return Kind::apolar;
  }
  if (name == "waterO") {
    return Kind::waterO;
  }
  if (name == "waterH") {
    return Kind::waterH;
  }
  if (name == "solvent") {
    return Kind::solvent;
  }
  throw std::invalid_argument("unknown site kind '" + std::string(name) + "'");
}

Family familyFromName(std::string_view name) {
  if (name == "waterIce") {
    return Family::waterIce;
  }
  if (name == "ionicLiquid") {
    return Family::ionicLiquid;
  }
  if (name == "moltenSalt") {
    return Family::moltenSalt;
  }
  if (name == "des") {
    return Family::des;
  }
  if (name == "electrolyte") {
    return Family::electrolyte;
  }
  if (name == "confinedIL") {
    return Family::confinedIL;
  }
  if (name == "confinedWater") {
    return Family::confinedWater;
  }
  if (name == "networkFormer") {
    return Family::networkFormer;
  }
  throw std::invalid_argument("unknown site family '" + std::string(name) +
                              "'");
}

bool isIonKind(Kind k) {
  return k == Kind::cationHead || k == Kind::anion;
}

bool matchesKind(Kind have, Kind want) {
  if (want == Kind::polar) {
    return have == Kind::cationHead || have == Kind::anion ||
           have == Kind::polar;
  }
  if (want == Kind::apolar) {
    return have == Kind::tail || have == Kind::apolar;
  }
  return have == want;
}

int indexOfAtom(const Cloud &src, int atomID) {
  const auto it = src.idIndexMap.find(atomID);
  if (it != src.idIndexMap.end()) {
    return it->second;
  }
  const int n = static_cast<int>(src.pts.size());
  for (int i = 0; i < n; ++i) {
    if (src.pts[static_cast<std::size_t>(i)].atomID == atomID) {
      return i;
    }
  }
  return -1;
}

} // namespace

Kind Table::of(const molSys::Point<double> &p) const {
  if (const auto it = atomOverride.find(p.atomID); it != atomOverride.end()) {
    return it->second;
  }
  return ofType(p.type);
}

Kind Table::ofType(int typeId) const {
  if (const auto it = typeToKind.find(typeId); it != typeToKind.end()) {
    return it->second;
  }
  return Kind::unspecified;
}

std::vector<int> indicesOf(const Cloud &yCloud, const Table &table,
                           Kind kind) {
  std::vector<int> out;
  const int n = static_cast<int>(yCloud.pts.size());
  out.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    if (matchesKind(table.of(yCloud.pts[static_cast<std::size_t>(i)]), kind)) {
      out.push_back(i);
    }
  }
  return out;
}

int lammpsTypeOfKind(const Table &table, Kind kind) {
  int found = 0;
  int typeId = -1;
  for (const auto &[id, mapped] : table.typeToKind) {
    if (mapped == kind) {
      ++found;
      typeId = id;
    }
  }
  if (found != 1) {
    throw std::runtime_error("site kind does not map to a unique LAMMPS type");
  }
  return typeId;
}

Cloud ionCloud(const Cloud &src, const Table &table) {
  Cloud out;
  out.box = src.box;
  out.boxLow = src.boxLow;
  out.currentFrame = src.currentFrame;

  auto &mutableSrc = const_cast<Cloud &>(src);
  const auto molMap = molSys::createMolIDAtomIDMultiMap(mutableSrc);

  std::vector<int> molOrder;
  std::unordered_set<int> seenMol;
  molOrder.reserve(src.pts.size());
  for (const auto &pt : src.pts) {
    if (seenMol.insert(pt.molID).second) {
      molOrder.push_back(pt.molID);
    }
  }

  for (const int molID : molOrder) {
    std::vector<int> ions;
    const auto range = molMap.equal_range(molID);
    for (auto it = range.first; it != range.second; ++it) {
      const int idx = indexOfAtom(src, it->second);
      if (idx < 0) {
        continue;
      }
      if (isIonKind(table.of(src.pts[static_cast<std::size_t>(idx)]))) {
        ions.push_back(idx);
      }
    }
    if (ions.empty()) {
      continue;
    }
    std::sort(ions.begin(), ions.end());
    const Kind ionKind = table.of(src.pts[static_cast<std::size_t>(ions.front())]);
    ions.erase(std::remove_if(ions.begin(), ions.end(),
                              [&](int idx) {
                                return table.of(src.pts[static_cast<std::size_t>(
                                           idx)]) != ionKind;
                              }),
               ions.end());
    if (ions.empty()) {
      continue;
    }

    molSys::Point<double> vertex = src.pts[static_cast<std::size_t>(ions.front())];
    vertex.type = (ionKind == Kind::cationHead) ? 1 : 2;
    vertex.molID = molID;
    if (ions.size() > 1) {
      const int ref = ions.front();
      double sx = src.pts[static_cast<std::size_t>(ref)].x;
      double sy = src.pts[static_cast<std::size_t>(ref)].y;
      double sz = src.pts[static_cast<std::size_t>(ref)].z;
      for (std::size_t k = 1; k < ions.size(); ++k) {
        const auto dr = gen::relDist(src, ions[k], ref);
        sx += src.pts[static_cast<std::size_t>(ref)].x + dr[0];
        sy += src.pts[static_cast<std::size_t>(ref)].y + dr[1];
        sz += src.pts[static_cast<std::size_t>(ref)].z + dr[2];
      }
      const double inv = 1.0 / static_cast<double>(ions.size());
      vertex.x = sx * inv;
      vertex.y = sy * inv;
      vertex.z = sz * inv;
    }

    const int outIdx = static_cast<int>(out.pts.size());
    out.idIndexMap[vertex.atomID] = outIdx;
    out.pts.push_back(std::move(vertex));
  }

  out.nop = static_cast<int>(out.pts.size());
  return out;
}

Table parseSiteSpec(std::string_view spec) {
  Table table;
  std::size_t start = 0;
  while (start <= spec.size()) {
    const std::size_t comma = spec.find(',', start);
    const auto raw =
        spec.substr(start, comma == std::string_view::npos
                               ? std::string_view::npos
                               : comma - start);
    start = (comma == std::string_view::npos) ? spec.size() + 1 : comma + 1;
    const std::string token = trim(raw);
    if (token.empty()) {
      continue;
    }
    const auto eq = token.find('=');
    if (eq == std::string::npos) {
      throw std::invalid_argument("site spec token '" + token +
                                  "' is not key=value");
    }
    const std::string key = trim(token.substr(0, eq));
    const std::string val = trim(token.substr(eq + 1));
    if (key.empty() || val.empty()) {
      throw std::invalid_argument("site spec token '" + token +
                                  "' is not key=value");
    }
    if (key == "family") {
      table.family = familyFromName(val);
      continue;
    }
    std::size_t consumed = 0;
    int typeId = 0;
    try {
      typeId = std::stoi(key, &consumed);
    } catch (const std::exception &) {
      throw std::invalid_argument("site spec type '" + key +
                                  "' is not an integer");
    }
    if (consumed != key.size()) {
      throw std::invalid_argument("site spec type '" + key +
                                  "' is not an integer");
    }
    table.typeToKind[typeId] = kindFromName(val);
  }
  return table;
}

IonEnvironment
ionEnvironment(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<bool> &iceFlag, const std::vector<int> &ionIndices,
               int waterType, double cutoff) {
  IonEnvironment out;
  const int n = yCloud.nop;
  std::vector<char> isIon(static_cast<std::size_t>(std::max(n, 0)), 0);
  for (int i : ionIndices) {
    if (i >= 0 && i < n) {
      isIon[static_cast<std::size_t>(i)] = 1;
    }
  }
  const double cut2 = cutoff * cutoff;
  for (int i : ionIndices) {
    if (i < 0 || i >= n) {
      continue;
    }
    int shell = 0;
    int labelled = 0;
    std::vector<int> members;
    for (int j = 0; j < n; j++) {
      if (j == i || isIon[static_cast<std::size_t>(j)]) {
        continue;
      }
      if (waterType != 0 && yCloud.pts[static_cast<std::size_t>(j)].type != waterType) {
        continue;
      }
      if (gen::periodicDistSq(yCloud, i, j) >= cut2) {
        continue;
      }
      ++shell;
      members.push_back(j);
      if (static_cast<std::size_t>(j) < iceFlag.size() && iceFlag[static_cast<std::size_t>(j)]) {
        ++labelled;
      }
    }
    IonState state = IonState::liquid;
    if (shell > 0 && labelled == shell) {
      state = IonState::ice;
      ++out.nIce;
    } else if (labelled > 0) {
      state = IonState::front;
      ++out.nFront;
    } else {
      ++out.nLiquid;
    }
    out.ion.push_back(i);
    out.shell.push_back(shell);
    out.iceFraction.push_back(shell > 0 ? static_cast<double>(labelled) / shell : 0.0);
    out.state.push_back(state);
    out.members.push_back(std::move(members));
  }
  return out;
}

std::vector<int> shellRingCensus(const std::vector<std::vector<int>> &rings,
                                 const std::vector<int> &shell, int maxRingSize) {
  std::vector<int> census(static_cast<std::size_t>(std::max(maxRingSize, 0)) + 1, 0);
  std::unordered_set<int> inShell(shell.begin(), shell.end());
  for (const auto &ring : rings) {
    if (ring.empty() || static_cast<int>(ring.size()) > maxRingSize) {
      continue;
    }
    for (int a : ring) {
      if (inShell.count(a)) {
        census[ring.size()] += 1;
        break;
      }
    }
  }
  return census;
}

ShellRings shellRings(const std::vector<std::vector<int>> &waterRings,
                      const std::vector<std::vector<int>> &nList, int ion,
                      const std::vector<int> &shell, int maxRingSize) {
  ShellRings out;
  out.census = shellRingCensus(waterRings, shell, maxRingSize);
  std::unordered_set<int> inShell(shell.begin(), shell.end());
  for (const auto &ring : waterRings) {
    if (ring.empty() || static_cast<int>(ring.size()) > maxRingSize) {
      continue;
    }
    bool allIn = true;
    for (int a : ring) {
      if (!inShell.count(a)) {
        allIn = false;
        break;
      }
    }
    if (allIn) {
      ++out.capped;
    }
  }
  if (ion < 0 || maxRingSize < 3 || shell.empty()) {
    return out;
  }
  int nVert = static_cast<int>(nList.size());
  for (int s : shell) {
    nVert = std::max(nVert, s + 1);
  }
  nVert = std::max(nVert, ion + 1);
  std::vector<std::vector<int>> aug(static_cast<std::size_t>(nVert));
  for (int i = 0; i < nVert; ++i) {
    aug[static_cast<std::size_t>(i)].push_back(i);
  }
  auto addEdge = [&](int a, int b) {
    if (a < 0 || b < 0 || a >= nVert || b >= nVert || a == b) {
      return;
    }
    auto &row = aug[static_cast<std::size_t>(a)];
    if (std::find(row.begin(), row.end(), b) == row.end()) {
      row.push_back(b);
    }
  };
  for (const auto &row : nList) {
    if (row.size() < 2) {
      continue;
    }
    const int i = row[0];
    for (std::size_t k = 1; k < row.size(); ++k) {
      addEdge(i, row[k]);
      addEdge(row[k], i);
    }
  }
  for (int s : shell) {
    addEdge(ion, s);
    addEdge(s, ion);
  }
  const auto plus = primitive::ringNetwork(aug, maxRingSize);
  for (const auto &ring : plus) {
    if (ring.empty() || static_cast<int>(ring.size()) > maxRingSize) {
      continue;
    }
    if (std::find(ring.begin(), ring.end(), ion) != ring.end()) {
      ++out.broken;
    }
  }
  return out;
}

namespace {
// minimum-image displacement r - p for a point r and a reference p
std::array<double, 3>
minImage(const molSys::PointCloud<molSys::Point<double>, double> &yCloud, double rx,
         double ry, double rz, double px, double py, double pz) {
  if (yCloud.box.size() >= 6) {
    return gen::triclinicMinImage(yCloud, rx, ry, rz, px, py, pz);
  }
  std::array<double, 3> dr = {rx - px, ry - py, rz - pz};
  for (int k = 0; k < 3; k++) {
    const double L = yCloud.box[static_cast<std::size_t>(k)];
    if (L > 0.0) {
      dr[static_cast<std::size_t>(k)] -= L * std::round(dr[static_cast<std::size_t>(k)] / L);
    }
  }
  return dr;
}
} // namespace

std::array<double, 3>
periodicCentroid(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                 const std::vector<int> &atoms) {
  std::array<double, 3> c = {0.0, 0.0, 0.0};
  if (atoms.empty()) {
    return c;
  }
  const auto &p0 = yCloud.pts[static_cast<std::size_t>(atoms.front())];
  for (int a : atoms) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(a)];
    const auto dr = minImage(yCloud, p.x, p.y, p.z, p0.x, p0.y, p0.z);
    for (int k = 0; k < 3; k++) {
      c[static_cast<std::size_t>(k)] += dr[static_cast<std::size_t>(k)];
    }
  }
  const double inv = 1.0 / static_cast<double>(atoms.size());
  c[0] = p0.x + c[0] * inv;
  c[1] = p0.y + c[1] * inv;
  c[2] = p0.z + c[2] * inv;
  return c;
}

GuestOccupancy
guestOccupancy(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
               const std::vector<std::vector<int>> &cages,
               const std::vector<int> &guestIndices, double radius) {
  GuestOccupancy out;
  out.guestsPerCage.assign(cages.size(), 0);
  std::vector<std::array<double, 3>> centres;
  centres.reserve(cages.size());
  for (const auto &cage : cages) {
    centres.push_back(periodicCentroid(yCloud, cage));
  }
  const double r2max = radius * radius;
  for (int g : guestIndices) {
    if (g < 0 || g >= yCloud.nop) {
      continue;
    }
    const auto &p = yCloud.pts[static_cast<std::size_t>(g)];
    int best = -1;
    double bestSq = r2max;
    for (std::size_t c = 0; c < centres.size(); c++) {
      const auto dr = minImage(yCloud, p.x, p.y, p.z, centres[c][0], centres[c][1], centres[c][2]);
      const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
      if (d2 <= bestSq) {
        bestSq = d2;
        best = static_cast<int>(c);
      }
    }
    out.cageOfGuest.push_back(best);
    out.centreDistance.push_back(best < 0 ? -1.0 : std::sqrt(bestSq));
    if (best < 0) {
      ++out.free;
    } else {
      ++out.guestsPerCage[static_cast<std::size_t>(best)];
    }
  }
  for (int n : out.guestsPerCage) {
    if (n > 0) {
      ++out.occupied;
    }
    if (n > 1) {
      ++out.multiply;
    }
  }
  return out;
}

namespace {

using Vec3 = std::array<double, 3>;

Vec3 vadd(const Vec3 &a, const Vec3 &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
Vec3 vsub(const Vec3 &a, const Vec3 &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
Vec3 vscale(const Vec3 &a, double s) {
  return {a[0] * s, a[1] * s, a[2] * s};
}
double vdot(const Vec3 &a, const Vec3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vec3 vcross(const Vec3 &a, const Vec3 &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}
double vnorm(const Vec3 &a) { return std::sqrt(vdot(a, a)); }

Vec3 unwrapAbout(const Cloud &yCloud, double x, double y, double z,
                 const Vec3 &c) {
  const auto dr = minImage(yCloud, x, y, z, c[0], c[1], c[2]);
  return vadd(c, dr);
}

struct Triangle {
  Vec3 a, b, c;
};

std::vector<Triangle> fanFace(const std::vector<Vec3> &verts) {
  std::vector<Triangle> out;
  if (verts.size() < 3) {
    return out;
  }
  Vec3 f = {0.0, 0.0, 0.0};
  for (const auto &v : verts) {
    f = vadd(f, v);
  }
  f = vscale(f, 1.0 / static_cast<double>(verts.size()));
  for (std::size_t i = 0; i < verts.size(); ++i) {
    out.push_back({f, verts[i], verts[(i + 1) % verts.size()]});
  }
  return out;
}

bool pointInTriangle(const Vec3 &p, const Triangle &t, double eps) {
  const Vec3 u = vsub(t.b, t.a);
  const Vec3 v = vsub(t.c, t.a);
  const Vec3 n = vcross(u, v);
  const double n2 = vdot(n, n);
  if (n2 < eps * eps) {
    return false;
  }
  if (std::fabs(vdot(vsub(p, t.a), n)) > eps * std::sqrt(n2)) {
    return false;
  }
  const Vec3 w = vsub(p, t.a);
  const double uu = vdot(u, u);
  const double uv = vdot(u, v);
  const double vv = vdot(v, v);
  const double wu = vdot(w, u);
  const double wv = vdot(w, v);
  const double det = uu * vv - uv * uv;
  if (std::fabs(det) < eps * eps) {
    return false;
  }
  const double s = (vv * wu - uv * wv) / det;
  const double r = (uu * wv - uv * wu) / det;
  return s >= -eps && r >= -eps && s + r <= 1.0 + eps;
}

enum class Hit { miss, edge, interior };

Hit rayTriangle(const Vec3 &orig, const Vec3 &dir, const Triangle &t,
                double eps) {
  const Vec3 e1 = vsub(t.b, t.a);
  const Vec3 e2 = vsub(t.c, t.a);
  const Vec3 pvec = vcross(dir, e2);
  const double det = vdot(e1, pvec);
  if (std::fabs(det) < eps) {
    return Hit::miss;
  }
  const double inv = 1.0 / det;
  const Vec3 tvec = vsub(orig, t.a);
  const double u = vdot(tvec, pvec) * inv;
  const Vec3 qvec = vcross(tvec, e1);
  const double v = vdot(dir, qvec) * inv;
  const double tt = vdot(e2, qvec) * inv;
  if (tt <= eps) {
    return Hit::miss;
  }
  if (u < -eps || v < -eps || u + v > 1.0 + eps) {
    return Hit::miss;
  }
  const bool onEdge = u <= eps || v <= eps || u + v >= 1.0 - eps;
  return onEdge ? Hit::edge : Hit::interior;
}

bool insideTriangles(const Vec3 &p, const std::vector<Triangle> &tris) {
  constexpr double eps = 1e-9;
  for (const auto &t : tris) {
    if (pointInTriangle(p, t, 1e-8)) {
      return true;
    }
  }
  for (int attempt = 0; attempt < 8; ++attempt) {
    Vec3 dir = {1.0, 1.0e-4 * static_cast<double>(attempt + 1),
                2.0e-4 * static_cast<double>(attempt + 1)};
    const double n = vnorm(dir);
    dir = vscale(dir, 1.0 / n);
    int hits = 0;
    bool grazed = false;
    for (const auto &t : tris) {
      const Hit h = rayTriangle(p, dir, t, eps);
      if (h == Hit::edge) {
        grazed = true;
        break;
      }
      if (h == Hit::interior) {
        ++hits;
      }
    }
    if (grazed) {
      continue;
    }
    return (hits % 2) == 1;
  }
  return false;
}

std::vector<Triangle>
cageTriangles(const Cloud &yCloud, const std::vector<int> &vertices,
              const std::vector<std::vector<int>> &faces, const Vec3 &centre) {
  std::vector<Triangle> tris;
  auto unwrapAtom = [&](int a) {
    const auto &p = yCloud.pts[static_cast<std::size_t>(a)];
    return unwrapAbout(yCloud, p.x, p.y, p.z, centre);
  };
  if (!faces.empty()) {
    for (const auto &face : faces) {
      std::vector<Vec3> verts;
      verts.reserve(face.size());
      for (int a : face) {
        if (a >= 0 && a < yCloud.nop) {
          verts.push_back(unwrapAtom(a));
        }
      }
      auto fan = fanFace(verts);
      tris.insert(tris.end(), fan.begin(), fan.end());
    }
    return tris;
  }
  // no faces: nothing to test
  (void)vertices;
  return tris;
}

void tallyOccupancy(GuestOccupancy &out) {
  out.occupied = 0;
  out.multiply = 0;
  out.free = 0;
  for (int n : out.guestsPerCage) {
    if (n > 0) {
      ++out.occupied;
    }
    if (n > 1) {
      ++out.multiply;
    }
  }
  for (int c : out.cageOfGuest) {
    if (c < 0) {
      ++out.free;
    }
  }
}

} // namespace

GuestOccupancy
guestOccupancyInside(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                     const std::vector<std::vector<int>> &cages,
                     const std::vector<std::vector<std::vector<int>>> &faces,
                     const std::vector<int> &guestIndices) {
  GuestOccupancy out;
  out.guestsPerCage.assign(cages.size(), 0);
  std::vector<Vec3> centres;
  std::vector<std::vector<Triangle>> meshes;
  centres.reserve(cages.size());
  meshes.reserve(cages.size());
  for (std::size_t c = 0; c < cages.size(); ++c) {
    centres.push_back(periodicCentroid(yCloud, cages[c]));
    const std::vector<std::vector<int>> empty;
    const auto &fc = c < faces.size() ? faces[c] : empty;
    meshes.push_back(cageTriangles(yCloud, cages[c], fc, centres.back()));
  }
  for (int g : guestIndices) {
    if (g < 0 || g >= yCloud.nop) {
      continue;
    }
    const auto &p = yCloud.pts[static_cast<std::size_t>(g)];
    int best = -1;
    double bestSq = 0.0;
    for (std::size_t c = 0; c < cages.size(); ++c) {
      if (meshes[c].empty()) {
        continue;
      }
      const Vec3 q = unwrapAbout(yCloud, p.x, p.y, p.z, centres[c]);
      if (!insideTriangles(q, meshes[c])) {
        continue;
      }
      const auto dr =
          minImage(yCloud, p.x, p.y, p.z, centres[c][0], centres[c][1],
                   centres[c][2]);
      const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
      if (best < 0 || d2 < bestSq) {
        best = static_cast<int>(c);
        bestSq = d2;
      }
    }
    out.cageOfGuest.push_back(best);
    out.centreDistance.push_back(best < 0 ? -1.0 : std::sqrt(bestSq));
    if (best >= 0) {
      ++out.guestsPerCage[static_cast<std::size_t>(best)];
    }
  }
  tallyOccupancy(out);
  return out;
}

DualOccupancy
guestOccupancyBoth(const molSys::PointCloud<molSys::Point<double>, double> &yCloud,
                   const std::vector<std::vector<int>> &cages,
                   const std::vector<std::vector<std::vector<int>>> &faces,
                   const std::vector<int> &guestIndices, double radius) {
  DualOccupancy out;
  out.radius = guestOccupancy(yCloud, cages, guestIndices, radius);
  out.inside = guestOccupancyInside(yCloud, cages, faces, guestIndices);
  return out;
}

} // namespace site
