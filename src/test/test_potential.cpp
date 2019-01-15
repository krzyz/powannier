#include <catch.hpp>
#include <vector>
#include "../powannier/builders/potentialbuilders.h"
#include "helpers.h"

TEST_CASE("Check Potential class", "[potential]") {
  const double pi = POWannier::pi;
  const double sqrt2 = std::sqrt(2);

  SECTION("Should not be able to create potential without basis") {
    auto builder = POWannier::Potential::create();
    REQUIRE_THROWS(builder.addFourierCoefficients());
    REQUIRE_THROWS(builder.evaluateFromFunction([](POWannier::Position x){return 0.0;}));
  }

  SECTION("Create simple 1d potentials") {
    std::vector<POWannier::Vector> basis({ {1} });

    auto&& func = [=](POWannier::Position x) {
      return std::cos(2 * pi * x[0]);
    };

    auto builder = POWannier::Potential::create();

    builder.setBasis(basis);

    std::unique_ptr<POWannier::Potential> V = nullptr;

    std::vector<POWannier::NPoint> indices = {{-1}, {1}};
    std::vector<POWannier::Complex> values = {0.5, 0.5};

    SECTION("by supplying all Fourier coefficients at once") {
      V = std::make_unique<POWannier::Potential>(builder
          .addFourierCoefficients()
          .provideMultiple(indices).values(values)
          .complete()
      );
    }

    SECTION("by supplying Fourier coefficients in parts") {
      auto x = builder.addFourierCoefficients();
      for (std::size_t i = 0; i < indices.size(); ++i) {
        x.provideOne(indices[i]).value(values[i]);
      }
      V = std::make_unique<POWannier::Potential>(builder
          .complete()
      );
    }

    SECTION("by supplying Fourier coefficients in parts") {
      V = std::make_unique<POWannier::Potential>(builder
          .evaluateFromFunction(func)
          .setCutoff(1)
          .setPrecision(1e-15)
          .requireRelativeErrorBelow(1e-14)
          .complete()
        );
    }

    for (auto x : std::vector<double>({-0.2, 0, 0.2})) {
      REQUIRE((*V)({x}) == Approx(func({x})).margin(1e-15));
    }
  }

  SECTION("Create inseparable 2D potential") {
    std::vector<POWannier::Vector> basis({ 
        {2, -2}, 
        {0, sqrt2} 
    });

    auto&& func = [=](POWannier::Position r) {
      auto x = r[0];
      auto y = r[1];
      return std::pow(std::sin(pi * (x + y) / sqrt2), 2) +
             std::pow(std::sin(0.5 * pi * x), 2) +
             std::pow(std::sin(pi * (x + 0.25)), 2);
    };

    auto V = POWannier::Potential::create()
          .setBasis(basis)
          .evaluateFromFunction(func)
          .setCutoff(2)
          .requireRelativeErrorBelow(1e-14)
          .complete();

    for (auto alpha1 : arma::linspace<arma::vec>(0, 1, 10)) {
      for (auto alpha2 : arma::linspace<arma::vec>(0, 1, 10)) {
        POWannier::Position r = alpha1 * basis[0] +
                                alpha2 * basis[1];
        REQUIRE(V(r) == Approx(func(r)).margin(1e-14));
      }
    }
  }
}
