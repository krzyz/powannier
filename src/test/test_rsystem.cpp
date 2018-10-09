#include <catch.hpp>

#include "../powannier/rsystem.h"
#include "../powannier/blochhelpers.h"

TEST_CASE("Check RSystem class", "[rsystem]") {
  const double pi = POWannier::pi;
  const double sqrt2 = std::sqrt(2);

  SECTION("Calculate rsystem for a given 1D potential") {
    double kl = pi;

    auto basis = std::vector<POWannier::Vector> { {pi/kl} };

    auto&& func = [=](POWannier::Position x) {
      return std::pow(std::sin(kl * x[0]), 2);
    };

    auto V = POWannier::Potential(basis, func, 1);

    auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 20, 8, 8);

    POWannier::RSystem rs(bs, {0});

    POWannier::NPoint n({0});
    auto wan = rs.getWannier(n);

    wan.drawB("wan1d", 20, POWannier::DrawRange::WholeLattice);
  }


  SECTION("Calculate rsystem for a given 2D potential") {
    double kl = pi;
    std::vector<POWannier::Vector> basis({ 
        {pi/kl, 0}, 
        {0, pi/kl} });

    auto&& func = [=](POWannier::Position r) {
      auto x = r[0];
      auto y = r[1];
      return std::pow(std::sin(kl * x), 2) * 
             std::pow(std::sin(kl * y), 2);
    };

    auto V = POWannier::Potential(basis, func, 1);

    auto N = 5;

    auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 10, N, 8);

    POWannier::RSystem rs(bs, {0});

    POWannier::NPoint n({0, 0});
    auto wan = rs.getWannier(n);

    //wan.drawB("wan2d", 100, POWannier::WholeLattice);
  }

  SECTION("Calculate rsystem for another 2D potential") {
    std::cout << "4 bands" << std::endl;
    double kl = pi;
    std::vector<POWannier::Vector> basis({ 
        {2*pi/kl, 0}, 
        {0, 2*pi/kl} });

    auto&& func = [=](POWannier::Position r) {
      auto x = r[0];
      auto y = r[1];
      return -std::pow(std::sin(kl * x), 2) * 
             std::pow(std::sin(kl * y), 2);
    };

    auto V = POWannier::Potential(basis, func, 2);
    REQUIRE(V.relative_error() == Approx(0).margin(1e-14));

    int N = 3;

    auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 10, N, 40);

    POWannier::RSystem rs(bs, {0,1,2,3}, {2, 2});

    POWannier::NPoint n({0, 0});
    auto wan = rs.getWannier(n, {0, 0});

    //wan.drawB("wan2d", 100, POWannier::WholeLattice);

/*
    for (int a = -1; a <= 1; ++a) {
      for (int b = -1; b <= 1; ++b) {
        POWannier::NPoint n({a, b});
        for (int ac = 0; ac < 2; ++ac) {
          for (int bc = 0; bc < 2; ++bc) {
            auto wan = rs.getWannier(n, {ac, bc});

            wan.drawB("wan2d" + std::to_string(a) + std::to_string(b)+std::to_string(ac) + std::to_string(bc), 100, POWannier::WholeLattice);
          }
        }
      }
    }
    */

 }

SECTION("Calculate rsystem for approximated 2D potential") {
    double kl = pi;
    std::vector<POWannier::Vector> basis({ 
        {pi/kl, 0}, 
        {0, pi/kl} });

    auto&& gauss1d = [=] (double x, double mu, double sigma) {
      return std::exp(-std::pow(x-mu,2) / (2*sigma*sigma) )/(std::sqrt(2*pi) * sigma);
    };

    auto&& gauss2d = [=] (POWannier::Position r, double mux, double muy, double sigma) {
      auto x = r[0];
      auto y = r[1];
      return gauss1d(x, mux, sigma) * gauss1d(y, muy, sigma);
    };


    auto&& func = [&](POWannier::Position r) {
      return -gauss2d(r, -0.25, -0.25, 0.1)
             -gauss2d(r, -0.25, 0.25, 0.1)
             -gauss2d(r, 0.25, 0, 0.1);
    };

    auto V = POWannier::Potential(basis, func, 2, 1e-12);

    int N = 3;

    auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 10, N, 40);

    POWannier::RSystem rs(bs, {0,1,2}, {3});

    POWannier::NPoint n({0, 0});
    auto wan1 = rs.getWannier(n, {0, 0});
    auto wan2 = rs.getWannier(n, {0, 1});
    auto wan3 = rs.getWannier(n, {0, 2});

    //drawBands(5, V, kl, 10, 1, "bands", 20, POWannier::DrawRange::BrillouinZone);

    wan1.drawB("wan2d3-1", 100, POWannier::DrawRange::WholeLattice);
    wan2.drawB("wan2d3-2", 100, POWannier::DrawRange::WholeLattice);
    wan3.drawB("wan2d3-3", 100, POWannier::DrawRange::WholeLattice);
 }

SECTION("Calculate rsystem for approximated 2D potential transposed") {
    double kl = pi;
    std::vector<POWannier::Vector> basis({ 
        {pi/kl, 0}, 
        {0, pi/kl} });

    auto&& gauss1d = [=] (double x, double mu, double sigma) {
      return std::exp(-std::pow(x-mu,2) / (2*sigma*sigma) )/(std::sqrt(2*pi) * sigma);
    };

    auto&& gauss2d = [=] (POWannier::Position r, double mux, double muy, double sigma) {
      auto x = r[1];
      auto y = r[0];
      return gauss1d(x, mux, sigma) * gauss1d(y, muy, sigma);
    };


    auto&& func = [&](POWannier::Position r) {
      return -gauss2d(r, -0.25, -0.25, 0.1)
             -gauss2d(r, -0.25, 0.25, 0.1)
             -gauss2d(r, 0.25, 0, 0.1);
    };

    auto V = POWannier::Potential(basis, func, 2, 1e-12);

    int N = 5;

    auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 10, N, 40);

    POWannier::RSystem rs(bs, {0,1,2}, {1, 1, 1});

    POWannier::NPoint n({0, 0});
    auto wan1 = rs.getWannier(n, {0, 0});
    auto wan2 = rs.getWannier(n, {1, 0});
    auto wan3 = rs.getWannier(n, {2, 0});

    wan1.drawB("wan2d3T-1", 100, POWannier::DrawRange::WholeLattice);
    wan2.drawB("wan2d3T-2", 100, POWannier::DrawRange::WholeLattice);
    wan3.drawB("wan2d3T-3", 100, POWannier::DrawRange::WholeLattice);
 }


  SECTION("Calculate rsystem for a given 3D potential") {
    double kl = pi;
    std::vector<POWannier::Vector> basis({ 
        {pi/kl, 0, 0}, 
        {0, pi/kl, 0},
        {0, 0, pi/kl} });

    auto&& func = [=](POWannier::Position r) {
      auto x = r[0];
      auto y = r[1];
      auto z = r[2];
      return std::pow(std::sin(kl * x), 2) * 
             std::pow(std::sin(kl * y), 2) *
             std::pow(std::sin(kl * z), 2);
    };

    auto V = POWannier::Potential(basis, func, 1);

    int N = 3;

    auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 1, N, 8);

    POWannier::RSystem rs(bs, {0});

    POWannier::NPoint n({0, 0, 0});
    auto wan = rs.getWannier(n);

    //wan.drawB("wan3d", 71, POWannier::WholeLattice);
 }

    //wan.drawB("wan3d", 71, POWannier::WholeLattice);
 

  /*
  SECTION("Calculate rsystem for a given 2D potential") {
    double kl = pi;
    std::vector<POWannier::Vector> basis({ 
        {2*pi/kl, -2*pi/kl}, 
        {0, sqrt2*pi/kl} });

    auto&& func = [=](POWannier::Position r) {
      auto x = r[0];
      auto y = r[1];
      return std::pow(std::sin(pi * (x + y) / sqrt2), 2) +
             std::pow(std::sin(0.5 * pi * x), 2) +
             std::pow(std::sin(pi * (x + 0.25)), 2);
    };

    auto V = POWannier::Potential(basis, func, 2);

    int N = 5;

    auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 4, N, 8);

    POWannier::RSystem rs(bs);

    POWannier::NPoint n({0, 0});

    auto wan = rs.getWannier(n);

    wan.drawB("wan2d", 20, POWannier::WholeLattice);
  }
    */
}
