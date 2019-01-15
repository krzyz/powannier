#include <catch.hpp>

#include "../powannier/rsystem.h"
#include "../powannier/blochhelpers.h"

TEST_CASE("Check RSystem class", "[rsystem]") {
  const double pi = POWannier::pi;

  SECTION("Calculate rsystem for a given 1D potential") {
    double kl = pi;

    auto basis = std::vector<POWannier::Vector> { {pi/kl} };

    auto&& func = [=](POWannier::Position x) {
      return std::pow(std::sin(kl * x[0]), 2);
    };

    auto V = POWannier::Potential::create()
      .setBasis(basis)
      .evaluateFromFunction(func)
      .setCutoff(1)
      .complete();

    auto bs = std::make_shared<POWannier::BlochSystem>(
      POWannier::BlochSystem::createWithPotential(V)
        .setLaserWavelength(kl)
        .setFourierCutoff(20)
        .setSystemSize(8)
        .setPotentialDepth(8)
        .complete()
    );

    POWannier::RSystem rs(bs, {0});

    POWannier::NPoint n({0});
    auto wan = rs.getWannier(n);
    //wan.draw("wan1d.dat", 100, POWannier::DrawRange::WholeLattice);
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

    auto V = POWannier::Potential::create()
      .setBasis(basis)
      .evaluateFromFunction(func)
      .setCutoff(1)
      .complete();

    auto N = 5;

    auto bs = std::make_shared<POWannier::BlochSystem>(
      POWannier::BlochSystem::createWithPotential(V)
        .setLaserWavelength(kl)
        .setFourierCutoff(10)
        .setSystemSize(N)
        .setPotentialDepth(8)
        .complete()
    );


    POWannier::RSystem rs(bs, {0});

    POWannier::NPoint n({0, 0});
    auto wan = rs.getWannier(n);
  }

  SECTION("Calculate rsystem for another 2D potential") {
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

    auto V = POWannier::Potential::create()
      .setBasis(basis)
      .evaluateFromFunction(func)
      .setCutoff(2)
      .complete();

    int N = 3;

    //auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 10, N, 40);
    auto bs = std::make_shared<POWannier::BlochSystem>(
      POWannier::BlochSystem::createWithPotential(V)
        .setLaserWavelength(kl)
        .setFourierCutoff(10)
        .setSystemSize(N)
        .setPotentialDepth(40)
        .complete()
    );


    POWannier::RSystem rs(bs, {0,1,2,3}, {2, 2});

    POWannier::NPoint n({0, 0});
    auto wan = rs.getWannier(n, {0, 0});
    /*
    auto eigs = rs.eigenvalues(n)
    for (auto x : eigs[0]) {
      std::cout << x << " " << std::endl;
    }
    for (auto x : eigs[1]) {
      std::cout << x << " " << std::endl;
    }
    */
    //wan.draw("wan2.dat", 50, POWannier::DrawRange::WholeLattice);


/*
    POWannier::RSystem rs(bs, {0,1,2,3}, {2, 2});

    POWannier::NPoint n({0, 0});
    auto wan = rs.getWannier(n, {0, 0});
    wan.draw("wan.dat", 50, POWannier::DrawRange::WholeLattice);
    */
 }

SECTION("Calculate rsystem for approximated 2D potential") {
    double kl = 1; 
    std::vector<POWannier::Vector> basis({ 
        {1/kl, 0}, 
        {0.5/kl, 0.5*std::sqrt(3)/kl} });

    auto&& gauss1d = [=] (double x, double mu, double sigma) {
      return std::exp(-std::pow(x-mu,2) / (2*sigma*sigma) )/(std::sqrt(2*pi) * sigma);
    };

    auto&& gauss2d = [=] (POWannier::Position r, double mux, double muy, double sigma) {
      auto x = r[0];
      auto y = r[1];
      return gauss1d(x, mux, sigma) * gauss1d(y, muy, sigma);
    };

    double a = 0.4;

    auto&& func = [&](POWannier::Position r) {
      r += POWannier::Position({0.5*a, 0.5*a/std::sqrt(3)});
      return -gauss2d(r, 0, 0, 0.1)
             -gauss2d(r, a, 0, 0.1)
             -gauss2d(r, 0.5*a, 0.5*std::sqrt(3)*a, 0.1);
    };

    auto V = POWannier::Potential::create()
      .setBasis(basis)
      .evaluateFromFunction(func)
      .setCutoff(3)
      .setPrecision(1e-12)
      .complete();


    int N = 3;

    //auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 5, N, 5);
    auto bs = std::make_shared<POWannier::BlochSystem>(
      POWannier::BlochSystem::createWithPotential(V)
        .setLaserWavelength(kl)
        .setFourierCutoff(5)
        .setSystemSize(N)
        .setPotentialDepth(5)
        .complete()
    );


    POWannier::RSystem rs(bs, {0,1,2}, {1, 1, 1});

    POWannier::NPoint n({0, 0});
    auto wan1 = rs.getWannier(n, {0, 0});
    auto wan2 = rs.getWannier(n, {1, 0});
    auto wan3 = rs.getWannier(n, {2, 0});
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

    auto V = POWannier::Potential::create()
      .setBasis(basis)
      .evaluateFromFunction(func)
      .setCutoff(1)
      .complete();


    int N = 3;

    //auto bs = std::make_shared<POWannier::BlochSystem>(V, kl, 1, N, 8);
    auto bs = std::make_shared<POWannier::BlochSystem>(
      POWannier::BlochSystem::createWithPotential(V)
        .setLaserWavelength(kl)
        .setFourierCutoff(2)
        .setSystemSize(N)
        .setPotentialDepth(8)
        .complete()
    );


    POWannier::RSystem rs(bs, {0});

    POWannier::NPoint n({0, 0, 0});
    auto wan = rs.getWannier(n);

    //wan.draw("wan3d.dat", 40, POWannier::DrawRange::WholeLattice);
 }
}
