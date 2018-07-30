#include <catch.hpp>
#include <vector>
#include "../powannier/potential.h"
#include "helpers.h"

TEST_CASE("Check Potential class", "[potential]") {
  const double pi = POWannier::pi;
  const double sqrt2 = std::sqrt(2);

  SECTION("Create simple 1D potential") {
    std::vector<POWannier::Vector> basis({ {1} });

    auto&& func = [=](POWannier::Position x) {
      return std::sin(2 * pi * x[0]);
    };

    POWannier::Potential V(basis, func, 1);
    for (auto x : std::vector<double>({-0.5, 0, 0.5})) {
      REQUIRE(V({x}) == Approx(func({x})).margin(1e-15));
    }
  }

  SECTION("Create inseparable 2D potential") {
    std::vector<POWannier::Vector> basis({ 
        {2, -2}, 
        {0, sqrt2} });

    auto&& func = [=](POWannier::Position r) {
      auto x = r[0];
      auto y = r[1];
      return std::pow(std::sin(pi * (x + y) / sqrt2), 2) +
             std::pow(std::sin(0.5 * pi * x), 2) +
             std::pow(std::sin(pi * (x + 0.25)), 2);
    };

    POWannier::Potential V(basis, func, 2);
    for (auto alpha1 : arma::linspace<arma::vec>(0, 1, 10)) {
      for (auto alpha2 : arma::linspace<arma::vec>(0, 1, 10)) {
        POWannier::Position r = alpha1 * basis[0] +
                                alpha2 * basis[1];
        REQUIRE(V({r}) == Approx(func({r})).margin(1e-14));
      }
    }
  }
}
