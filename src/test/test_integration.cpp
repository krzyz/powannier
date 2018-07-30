#include <catch.hpp>
#include <complex>
#include <cmath>
#include <functional>
#include <vector>
#include <armadillo>
#include "../util/cubatureintegration.h"
#include "helpers.h"

namespace CubatureIntegrationTest {
  double tangent(std::vector<double> x) {
    return std::tan(x[0]);
  }

  struct TestFunctor {
    double operator()(std::vector<double> x) {
      return std::log(x[0]);
    }
  };

  template <class IntegrationProvider>
  class IntegrationTestClass {
    public:
      IntegrationTestClass() {}
      void testMethod() {
        SECTION("Check if it implements IntegrationInterface") {
          REQUIRE(std::is_base_of<POWannier::IntegrationInterface, IntegrationProvider>::value == true);
        }

        SECTION("Calculate simple 1D integral using lambda expression") {
          std::vector<double> xmin = {0};
          std::vector<double> xmax = {1};
          double result = IntegrationProvider::integrate(
              [] (std::vector<double> x) {
                  auto value = std::sin(2*M_PI*x[0]) * std::sin(2*M_PI*x[0]);
                  return value;
              }, xmin, xmax);

          REQUIRE(result == Approx(0.5));
        }

        SECTION("Calculate simple 2d integral using lambda expression") {
          std::vector<double> xmin = {0, 0};
          std::vector<double> xmax = {1, 1};
          double result = IntegrationProvider::integrate(
              [] (std::vector<double> x) {
                  auto value = std::pow(std::sin(2*M_PI*x[0]), 2) *
                      std::pow(std::sin(2*M_PI*x[1]), 2);
                  return value;
               }, xmin, xmax);

          REQUIRE(result == Approx(0.25));
        }


        SECTION("Calculate integral using a function in a global namespace") {
          std::vector<double> xmin = {0};
          std::vector<double> xmax = {M_PI/4};
          double result = IntegrationProvider::integrate(
              &CubatureIntegrationTest::tangent, xmin, xmax);

          REQUIRE(result == Approx(0.5 * std::log(2)));
        }

        SECTION("Calculate integral using a functor") {
          std::vector<double> xmin = {1};
          std::vector<double> xmax = {std::exp(1.0)};
          CubatureIntegrationTest::TestFunctor testfunctor;
          double result = IntegrationProvider::integrate(
              std::ref(testfunctor), xmin, xmax);

          REQUIRE(result == Approx(1));
        }

        SECTION("Calculate two-dimensional integral") {
          std::vector<double> xmin = {0, 0};
          std::vector<double> xmax = {1, 1};
          double result = IntegrationProvider::integrate(
              [] (std::vector<double> x) {
                  auto value = std::sin(M_PI * (x[0] + x[1])) *
                    std::sin(M_PI * (x[0] - x[1])) *
                    std::cos(2 * M_PI * x[0]);
                  return value;
              }, xmin, xmax);
          REQUIRE(result == Approx(-0.25));
        }

     SECTION("Calculate complex integral") {
          std::vector<double> xmin = {0};
          std::vector<double> xmax = {1};
          std::complex<double> result = IntegrationProvider::integrate(
              [] (std::vector<double> x) {
                  auto value = std::sin(2 * M_PI * x[0]) *
                  std::exp(std::complex<double>(0, 2.0 * M_PI * x[0]));
                  return value;
              }, xmin, xmax);
          REQUIRE(std::real(result) == Approx(0).margin(1e-15));
          REQUIRE(std::imag(result) == Approx(0.5));
        }

     SECTION("Calculate complex two-dimensional integral") {
          arma::rowvec xmin(2, arma::fill::zeros);
          arma::rowvec xmax(2, arma::fill::ones);
          std::complex<double> result = IntegrationProvider::integrate(
              [=] (arma::rowvec r) {
                  double x = 2 * r[0];
                  double y = -2 * r[0] + std::sqrt(2) * r[1];
                  double an = -2 * r[0];
                  auto value = (std::pow(std::sin(M_PI * (x + y) / std::sqrt(2)), 2) +
                        std::pow(std::sin(0.5 * M_PI * x), 2) +
                        std::pow(std::sin(M_PI * (x + 0.25)), 2) )* 
                        std::exp(std::complex<double>(0, 2 * M_PI * an ));
                  return value;
              }, xmin, xmax);
          REQUIRE(std::real(result) == Approx(0).margin(1e-15));
          REQUIRE(std::imag(result) == Approx(-0.25));
        }
      }

  };

}


METHOD_AS_TEST_CASE(CubatureIntegrationTest::IntegrationTestClass<POWannier::CubatureIntegration>::testMethod,
    "Calculating integrals using cubature", "[integral, cubature]");
