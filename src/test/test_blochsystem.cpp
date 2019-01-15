#include <catch.hpp>
#include <vector>
#include <cmakedefines.h>
#include "../powannier/builders/potentialbuilders.h"
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

    auto V = POWannier::Potential::create()
      .setBasis(basis)
      .evaluateFromFunction(func)
      .setCutoff(1)
      .complete();

    auto bs = POWannier::BlochSystem::createWithPotential(V)
      .setLaserWavelength(kl)
      .setFourierCutoff(5)
      .setSystemSize(20)
      .setPotentialDepth(8)
      .complete();

    auto ms = std::vector<POWannier::NPoint> {{0}, {9}};
    auto bands = std::vector<int> {0, 1};
    auto i = 0;
    for (auto m : ms) {

      arma::cx_vec bf(5);
      auto x = arma::linspace<arma::vec>(0, 10, bf.n_elem);

      for (std::size_t j = 0; j < bf.n_elem; ++j) {
        bf(j) = bs.bloch(m, {x[j]}, bands[i]);
      }

      arma::cx_vec bfSaved;
      bfSaved.load(CMakeDefines::test_data_directory + "/blochFunction1d_" + std::to_string(i), arma::raw_ascii);

      for (std::size_t j = 0; j < bfSaved.n_elem; ++j) {
        REQUIRE(bf(j).real() == Approx(bfSaved(j).real()).margin(1e-14));
        REQUIRE(bf(j).imag() == Approx(bfSaved(j).imag()).margin(1e-14));
      }

      //bf.save("blochFunction1d_" + std::to_string(i), arma::raw_ascii);
      ++i;
    }
  }
}
