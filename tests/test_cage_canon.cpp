#include <catch2/catch_test_macros.hpp>

#include <cage_canon.hpp>

#include <thread>
#include <vector>

TEST_CASE("nautyAvailable matches the compile-time flag", "[nauty]") {
#ifdef SEAMS_HAS_NAUTY
  REQUIRE(cage::nautyAvailable());
#else
  REQUIRE_FALSE(cage::nautyAvailable());
  REQUIRE(cage::canonicalCertificate({{0, 1, 2}}).empty());
  REQUIRE(cage::canonicalCertificate({}).empty());
  REQUIRE_FALSE(cage::isHexagonalPrism({{0, 1, 2, 3, 4, 5}}));
  REQUIRE_FALSE(cage::isHexagonalPrism({}));
  REQUIRE_FALSE(cage::sameCertificate({{0, 1, 2}}, {{0, 1, 2}}));
  REQUIRE(cage::canonicalCertificateRooted(0, {}, -1).empty());
  REQUIRE(cage::canonicalCertificateRooted(2, {{0, 0}, {-1, 0}, {0, 3}}, 0).empty());
  REQUIRE(cage::canonicalCertificateColoured(1, {{0, 1}}, {0}, -1).empty());
#endif
}

#ifdef SEAMS_HAS_NAUTY

static std::vector<std::vector<int>> hexPrismRings() {
  return {{0, 1, 2, 3, 4, 5},
          {6, 7, 8, 9, 10, 11},
          {0, 1, 7, 6},
          {1, 2, 8, 7},
          {2, 3, 9, 8},
          {3, 4, 10, 9},
          {4, 5, 11, 10},
          {5, 0, 6, 11}};
}

static std::vector<std::vector<int>> relabel(const std::vector<std::vector<int>> &rings,
                                             int n) {
  std::vector<std::vector<int>> out = rings;
  for (auto &ring : out) {
    for (int &a : ring) {
      a = n - 1 - a;
    }
  }
  return out;
}

TEST_CASE("hexagonal prism certificate is stable under relabel", "[nauty]") {
  const auto rings = hexPrismRings();
  REQUIRE(cage::isHexagonalPrism(rings));
  REQUIRE(cage::sameCertificate(rings, relabel(rings, 12)));
  REQUIRE_FALSE(cage::isHexagonalPrism({{0, 1, 2, 3, 4, 5}}));
}

TEST_CASE("an empty ring set has an empty certificate", "[nauty]") {
  REQUIRE(cage::canonicalCertificate({}).empty());
  REQUIRE_FALSE(cage::isHexagonalPrism({}));
  REQUIRE_FALSE(cage::sameCertificate({}, {}));
  REQUIRE_FALSE(cage::sameCertificate({}, hexPrismRings()));
}

TEST_CASE("degenerate rings keep the underlying edge", "[nauty]") {
  const auto lone = cage::canonicalCertificate({{7}});
  REQUIRE_FALSE(lone.empty());
  REQUIRE_FALSE(cage::isHexagonalPrism({{7}}));
  const auto edge = cage::canonicalCertificate({{0, 1}});
  REQUIRE(cage::canonicalCertificate({{0, 0, 1}, {1, 0}}) == edge);
  REQUIRE(cage::canonicalCertificate({{0, 1}, {0}}) == edge);
  REQUIRE(cage::canonicalCertificate({{0, 1, 2}, {0, 1, 2}}) ==
          cage::canonicalCertificate({{0, 1, 2}}));
  REQUIRE_FALSE(cage::sameCertificate({{0, 1}}, {{0, 1, 2}}));
  REQUIRE(cage::canonicalCertificate({{0, 1, 2}}) !=
          cage::canonicalCertificate({{0}, {1}, {2}}));
}

TEST_CASE("rooted and coloured certificates follow the partition", "[nauty]") {
  const std::vector<std::pair<int, int>> path = {{0, 1}, {1, 2}};
  const auto end = cage::canonicalCertificateRooted(3, path, 0);
  const auto mid = cage::canonicalCertificateRooted(3, path, 1);
  const auto plain = cage::canonicalCertificateRooted(3, path, -1);
  REQUIRE_FALSE(end.empty());
  REQUIRE(end != mid);
  REQUIRE(cage::canonicalCertificateRooted(3, path, 3) == plain);
  REQUIRE(cage::canonicalCertificateRooted(0, path, 0).empty());
  REQUIRE(cage::canonicalCertificateRooted(-1, path, 0).empty());
  REQUIRE(cage::canonicalCertificateRooted(65, path, 0).empty());

  const std::vector<std::pair<int, int>> messy = {
      {-1, 0}, {0, -1}, {3, 0}, {0, 3}, {1, 1}, {0, 1}, {1, 0}, {0, 1}};
  REQUIRE(cage::canonicalCertificateRooted(3, messy, -1) ==
          cage::canonicalCertificateRooted(3, {{0, 1}}, -1));

  const auto ends = cage::canonicalCertificateColoured(3, path, {1, 0, 1}, -1);
  const auto midColour =
      cage::canonicalCertificateColoured(3, path, {0, 1, 0}, -1);
  REQUIRE_FALSE(ends.empty());
  REQUIRE(ends != midColour);
  const auto rootedColour =
      cage::canonicalCertificateColoured(3, path, {0, 0, 1}, 0);
  REQUIRE_FALSE(rootedColour.empty());
  REQUIRE(rootedColour != ends);
  REQUIRE(cage::canonicalCertificateColoured(3, messy, {0, 0, 1}, 0) ==
          cage::canonicalCertificateColoured(3, {{0, 1}}, {0, 0, 1}, 0));
  // A root outside the vertex range is not its own cell.
  REQUIRE(cage::canonicalCertificateColoured(3, path, {0, 0, 1}, 3) ==
          cage::canonicalCertificateColoured(3, path, {0, 0, 1}, -1));
  REQUIRE(cage::canonicalCertificateColoured(65, path, {0, 1, 2}, 0).empty());
  REQUIRE(cage::canonicalCertificateColoured(0, {}, {}, -1).empty());
}

TEST_CASE("two workers agree on one certificate", "[nauty]") {
  const auto rings = hexPrismRings();
  const std::string once = cage::canonicalCertificate(rings);
  REQUIRE_FALSE(once.empty());
  std::string left;
  std::string right;
  std::thread a([&] {
    for (int i = 0; i < 40; i++) {
      left = cage::canonicalCertificate(rings);
    }
  });
  std::thread b([&] {
    for (int i = 0; i < 40; i++) {
      right = cage::canonicalCertificateRooted(3, {{0, 1}, {1, 2}}, 0);
    }
  });
  a.join();
  b.join();
  REQUIRE(left == once);
  REQUIRE(right == cage::canonicalCertificateRooted(3, {{0, 1}, {1, 2}}, 0));
  REQUIRE(cage::isHexagonalPrism(rings));
}

TEST_CASE("no self-edge reaches the certificate", "[nauty]") {
  REQUIRE(cage::canonicalCertificate({{0, 0, 1}}) ==
          cage::canonicalCertificate({{0, 1}}));
  auto rings = hexPrismRings();
  rings.push_back({0, 0});
  REQUIRE(cage::isHexagonalPrism(rings));
  const std::vector<std::pair<int, int>> clean = {{0, 1}};
  const std::vector<std::pair<int, int>> withSelf = {{0, 0}, {1, 1}, {0, 1}};
  const std::vector<std::pair<int, int>> onlySelf = {{0, 0}, {1, 1}};
  REQUIRE(cage::canonicalCertificateRooted(2, withSelf, -1) ==
          cage::canonicalCertificateRooted(2, clean, -1));
  REQUIRE(cage::canonicalCertificateRooted(2, onlySelf, -1) ==
          cage::canonicalCertificateRooted(2, {}, -1));
  REQUIRE(cage::canonicalCertificateColoured(2, withSelf, {0, 1}, 0) ==
          cage::canonicalCertificateColoured(2, clean, {0, 1}, 0));
  REQUIRE(cage::canonicalCertificateColoured(2, onlySelf, {0, 1}, -1) ==
          cage::canonicalCertificateColoured(2, {}, {0, 1}, -1));
}

#endif
