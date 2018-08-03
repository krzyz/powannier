#include <catch.hpp>
#include <vector>
#include <iomanip>
#include "../powannier/potential.h"
#include "../powannier/blochsystem.h"
#include "../powannier/aliases.h"

TEST_CASE("Check BlochSystem class") {
  const double pi = POWannier::pi;

  SECTION("Calculate bloch function for a given 1D potential") {
    double kl = pi;

    std::vector<POWannier::Vector> basis({ {pi/kl} });

    auto&& func = [=](POWannier::Position x) {
      return std::pow(std::sin(kl * x[0]), 2);
    };

    auto V = std::make_shared<POWannier::Potential>(basis, func, 1);

    POWannier::BlochSystem bs(V, kl, 5, 20, 8);
    std::vector<POWannier::NPoint> ms = {{0}, {9}};
    std::vector<int> bands = {0, 1};
    int i = 0;
    for (auto m : ms) {

      arma::cx_vec bf(5);
      auto x = arma::linspace<arma::vec>(0, 10, bf.n_elem);

      for (int j = 0; j < bf.n_elem; ++j) {
        bf(j) = bs.bloch(m, {x[j]}, bands[i]);
      }

      arma::cx_vec bfSaved;
      bfSaved.load("blochFunction1d_" + std::to_string(i), arma::raw_ascii);

      for (int i = 0; i < bfSaved.n_elem; ++i) {
        REQUIRE(bf(i).real() == Approx(bfSaved(i).real()).margin(1e-14));
        REQUIRE(bf(i).imag() == Approx(bfSaved(i).imag()).margin(1e-14));
      }

      //bf.save("blochFunction1d_" + std::to_string(i), arma::raw_ascii);
      ++i;
    }
  }
}