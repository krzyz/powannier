#include <catch.hpp>

#include "../powannier/rsystem.h"

TEST_CASE("Check RSystem class") {
  const double pi = POWannier::pi;

  SECTION("Calculate rsystem for a given 1D potential") {
    double kl = pi;

    std::vector<POWannier::Vector> basis({ {pi/kl} });

    auto&& func = [=](POWannier::Position x) {
      return std::pow(std::sin(kl * x[0]), 2);
    };

    auto V = std::make_shared<POWannier::Potential>(basis, func, 1);

    auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 20, 8, 8);

    POWannier::RSystem rs(bs);

    POWannier::NPoint n({0});
    auto wan = rs.getWannier(n);

    arma::vec xs = arma::linspace<arma::vec>(-3, 3);
    for (auto r : xs) {
      std::cout << r << "\t" << wan({r}) << std::endl;
    }
  }
}
