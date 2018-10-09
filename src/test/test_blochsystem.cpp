#include <catch.hpp>
#include <vector>
#include "../powannier/potential.h"
#include "../powannier/blochsystem.h"
#include "../powannier/aliases.h"

TEST_CASE("Check BlochSystem class", "[bloch]") {
  const auto pi = POWannier::pi;

  SECTION("Calculate bloch function for a given 1D potential") {
    auto kl = pi;

    std::vector<POWannier::Vector> basis({ {pi/kl} });

    auto&& func = [=](POWannier::Position x) {
      return std::pow(std::sin(kl * x[0]), 2);
    };

    auto V = POWannier::Potential(basis, func, 1);

    auto bs = POWannier::BlochSystem (V, kl, 5, 20, 8);
    auto ms = std::vector<POWannier::NPoint> {{0}, {9}};
    auto bands = std::vector<int> {0, 1};
    auto i = 0;
    for (auto m : ms) {

      arma::cx_vec bf(5);
      auto x = arma::linspace<arma::vec>(0, 10, bf.n_elem);

      for (auto j = 0; j < bf.n_elem; ++j) {
        bf(j) = bs.bloch(m, {x[j]}, bands[i]);
      }

      arma::cx_vec bfSaved;
      bfSaved.load("blochFunction1d_" + std::to_string(i), arma::raw_ascii);

      for (auto i = 0; i < bfSaved.n_elem; ++i) {
        REQUIRE(bf(i).real() == Approx(bfSaved(i).real()).margin(1e-14));
        REQUIRE(bf(i).imag() == Approx(bfSaved(i).imag()).margin(1e-14));
      }

      //bf.save("blochFunction1d_" + std::to_string(i), arma::raw_ascii);
      ++i;
    }
  }
}
