#include <catch.hpp>
#include <vector>
#include "../powannier/potential.h"

TEST_CASE("Check Potential class", "[potential]") {
  SECTION("Creating simple 1D potential") {
    std::vector<POWannier::Vector> basis({1});
    auto&& func = [](POWannier::Position x) {
      return sin(x[0]);
    };
    POWannier::Potential V(basis, func, 1);
    for (auto x : std::vector<double>({-0.5, 0, 0.5})) {
      REQUIRE(V({x}) == func({x}));
    }
  }
}