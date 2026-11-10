#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sphericart_ylm.hpp>
#include <steinhardt_device.hpp>

#include <cmath>
#include <vector>

TEST_CASE("sphericart_ylm::available matches the compile-time flag",
          "[sphericart]") {
#ifdef SEAMS_HAS_SPHERICART
  REQUIRE(seams::sphericart_ylm::available());
#else
  REQUIRE_FALSE(seams::sphericart_ylm::available());
  double dummy = 0.0;
  REQUIRE(seams::sphericart_ylm::ylmCartesian(6, &dummy, 0, &dummy) == 0);
  const double xyz[3] = {0.3, 0.4, 0.8660254037844386};
  double out = 0.0;
  REQUIRE(seams::sphericart_ylm::ylmCartesian(6, xyz, 1, &out) == 1);
  REQUIRE(seams::sphericart_ylm::ylmCartesian(5, xyz, 1, &out) == 1);
#endif
}

#ifdef SEAMS_HAS_SPHERICART

TEST_CASE("sphericart_ylm accepts l=12", "[sphericart]") {
  const double xyz[3] = {0.3, 0.4, 0.8660254037844386};
  double cart[50];
  REQUIRE(seams::sphericart_ylm::ylmCartesian(12, xyz, 1, cart) == 0);
  double mod = 0.0;
  for (int m = 0; m < 25; m++) {
    mod += cart[2 * m] * cart[2 * m] + cart[2 * m + 1] * cart[2 * m + 1];
  }
  REQUIRE(mod > 0.0);
}

TEST_CASE("sphericart real Ylm have the same l=6 modulus as ylmAll",
          "[sphericart]") {
  const double xyz[3] = {0.3, 0.4, 0.8660254037844386};
  double cart[26];
  REQUIRE(seams::sphericart_ylm::ylmCartesian(6, xyz, 1, cart) == 0);
  const double r = std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
  const double phi = std::atan2(xyz[0], xyz[1]);
  const double theta = std::acos(xyz[2] / r);
  double closed[26];
  seams::steinhardt::ylmAll(6, theta, phi, closed);
  double mCart = 0.0;
  double mClosed = 0.0;
  for (int m = 0; m < 13; m++) {
    mCart += cart[2 * m] * cart[2 * m] + cart[2 * m + 1] * cart[2 * m + 1];
    mClosed +=
        closed[2 * m] * closed[2 * m] + closed[2 * m + 1] * closed[2 * m + 1];
  }
  REQUIRE_THAT(mCart, Catch::Matchers::WithinRel(mClosed, 1e-6));
}

TEST_CASE("sphericart accepts the implemented orders and rejects the rest",
          "[sphericart]") {
  REQUIRE(seams::sphericart_ylm::ylmCartesian(6, nullptr, 0, nullptr) == 0);
  const double xyz[3] = {0.3, 0.4, 0.8660254037844386};
  double scratch = 0.0;
  REQUIRE(seams::sphericart_ylm::ylmCartesian(5, xyz, 1, &scratch) == 1);
  REQUIRE(seams::sphericart_ylm::ylmCartesian(6, nullptr, -1, nullptr) == 1);
  for (int order : {3, 4, 8}) {
    std::vector<double> cart(static_cast<size_t>(2 * order + 1) * 2, -1.0);
    REQUIRE(seams::sphericart_ylm::ylmCartesian(order, xyz, 1, cart.data()) == 0);
    double mod = 0.0;
    for (size_t i = 0; i < cart.size(); i += 2) {
      REQUIRE(cart[i + 1] == 0.0);
      mod += cart[i] * cart[i];
    }
    REQUIRE(mod > 0.0);
  }
  const double pair[6] = {0.3, 0.4, 0.8660254037844386, 0.0, 0.0, 1.0};
  std::vector<double> both(26 * 2, -1.0);
  REQUIRE(seams::sphericart_ylm::ylmCartesian(6, pair, 2, both.data()) == 0);
  for (int row = 0; row < 2; row++) {
    double mod = 0.0;
    for (int m = 0; m < 13; m++) {
      const double re = both[static_cast<size_t>(row * 26 + 2 * m)];
      const double im = both[static_cast<size_t>(row * 26 + 2 * m + 1)];
      REQUIRE(im == 0.0);
      mod += re * re;
    }
    REQUIRE(mod > 0.0);
  }
}

#endif
