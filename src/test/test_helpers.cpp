#include <catch.hpp>
#include "../powannier/helpers.h"

TEST_CASE("Chop function works correctly", "[helpers, chop]") {

  std::vector<double> numbers = {-3, -0.5, -1e-15, 0, 1e-16, 0.5, 3};
  std::vector<double> numbersChopped = {-3, -0.5, 0, 0, 0, 0.5, 3};

  SECTION("Chopping real numbers") {
    for (size_t i = 0; i < numbers.size(); ++i) {
      REQUIRE(POWannier::chop(numbers[i], 1e-14) == numbersChopped[i]);
    }
  }

  SECTION("Chopping complex numbers") {
    std::vector<std::complex<double>> numbersComplex;
    std::vector<std::complex<double>> numbersComplexChopped;
    for (size_t i = 0; i < numbers.size(); ++i) {
      for (size_t j = 0; j < numbers.size(); ++j) {
        numbersComplex.push_back({numbers[i], numbers[j]});
        numbersComplexChopped.push_back({numbersChopped[i], numbersChopped[j]});
      }
    }

    for (size_t i = 0; i < numbersComplex.size(); ++i) {
      REQUIRE(POWannier::chop(numbersComplex[i], 1e-14) == numbersComplexChopped[i]);
    }
  }
}

TEST_CASE("Check generation of n-dimensional grid of uniformly distributed points with given cutoff", "[helpers, nspace]") {
  SECTION("Create 2-dimensional point grid") {
    std::vector<POWannier::NPoint> npoints = {
      {-1, -1}, {-1, 0}, {-1, 1},
      {0, -1}, {0, 0}, {0, 1},
      {1, -1}, {1, 0}, {1, 1}};
    int cutoff = 1;
    int dim = 2;
    std::vector<POWannier::NPoint> npointsGen = POWannier::nspace(cutoff, dim);
    REQUIRE(npoints.size() == npointsGen.size());
    for (size_t i = 0; i < npoints.size(); ++i) {
      REQUIRE(arma::all(npoints[i] == npointsGen[i]));
      REQUIRE(POWannier::nIndex(npoints[i], cutoff) == i);
    }
  }
}
