#include <catch.hpp>
#include "../powannier/helpers.h"

TEST_CASE("Chop function works correctly", "[helpers, chop]") {
  SECTION("Create simple 1D potential") {
    std::vector<POWannier::NPoint> npoints = {
      {-1, -1}, {-1, 0}, {-1, 1},
      {0, -1}, {0, 0}, {0, 1},
      {1, -1}, {1, 0}, {1, 1}};
    std::vector<POWannier::NPoint> npointsGen = POWannier::nspace(1, 2);
    REQUIRE(npoints.size() == npointsGen.size());
    for (int i = 0; i < npoints.size(); ++i) {
      REQUIRE(arma::all(npoints[i] == npointsGen[i]));
    }
  }
}

TEST_CASE("Check generation of n-dimensional uniformly distributed points with given cutoff", "[helpers, nspace]") {
  SECTION("Create simple 1D potential") {
    std::vector<POWannier::NPoint> npoints = {
      {-1, -1}, {-1, 0}, {-1, 1},
      {0, -1}, {0, 0}, {0, 1},
      {1, -1}, {1, 0}, {1, 1}};
    std::vector<POWannier::NPoint> npointsGen = POWannier::nspace(1, 2);
    REQUIRE(npoints.size() == npointsGen.size());
    for (int i = 0; i < npoints.size(); ++i) {
      REQUIRE(arma::all(npoints[i] == npointsGen[i]));
    }
  }
}
