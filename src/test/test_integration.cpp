#include <catch.hpp>
#include <complex>
#include <cmath>
#include <functional>
#include <vector>
#include "../util/cubatureintegration.h"

namespace CubatureIntegrationTest {
  double tangent(std::vector<double> x) {
    return std::tan(x[0]);
  }

  struct TestFunctor {
    double operator()(std::vector<double> x) {
      return std::log(x[0]);
    }
  };
}

template <class IntegrationProvider>
class IntegrationTestClass {
  public:
    IntegrationTestClass() {}
    void testMethod() {

      SECTION("Check if it implements IntegrationInterface") {
        REQUIRE(std::is_base_of<POWannier::IntegrationInterface, IntegrationProvider>::value == true);
      }
      SECTION("Calculate integral using lambda expression") {
        std::vector<double> xmin = {0};
        std::vector<double> xmax = {1};
        double result = IntegrationProvider::integrate(
            [] (std::vector<double> x) {
                return std::sin(2*M_PI*x[0]) * std::sin(2*M_PI*x[0]);
            }, xmin, xmax);

        REQUIRE(result == 0.5);
      }

      SECTION("Calculate integral using a function in a global namespace") {
        std::vector<double> xmin = {0};
        std::vector<double> xmax = {M_PI/4};
        double result = IntegrationProvider::integrate(
            &CubatureIntegrationTest::tangent, xmin, xmax);

        REQUIRE(result == 0.5 * std::log(2));
      }

      SECTION("Calculate integral using a functor") {
        std::vector<double> xmin = {1};
        std::vector<double> xmax = {std::exp(1.0)};
        CubatureIntegrationTest::TestFunctor testfunctor;
        double result = IntegrationProvider::integrate(
            std::ref(testfunctor), xmin, xmax);

        REQUIRE(result == 1);
      }

      SECTION("Calculate two-dimensional integral") {
        std::vector<double> xmin = {0, 0};
        std::vector<double> xmax = {1, 1};
        double result = IntegrationProvider::integrate(
            [] (std::vector<double> x) {
                return std::sin(M_PI * (x[0] + x[1])) *
                  std::sin(M_PI * (x[0] - x[1])) *
                  std::cos(2 * M_PI * x[0]);
            }, xmin, xmax);
        REQUIRE(result == -0.25);
      }

   SECTION("Calculate complex integral") {
          std::vector<double> xmin = {0};
          std::vector<double> xmax = {1};
          std::complex<double> result = IntegrationProvider::integrate(
              [] (std::vector<double> x) {
                  return std::sin(2 * M_PI * x[0]) *
                  std::exp(std::complex<double>(0, 2.0 * M_PI * x[0]));
              }, xmin, xmax);
          REQUIRE(result == std::complex<double>(0, 0.5));
        }
    }
};

METHOD_AS_TEST_CASE(IntegrationTestClass<POWannier::CubatureIntegration>::testMethod,
    "Calculating integrals using cubature", "[integral, cubature]");
