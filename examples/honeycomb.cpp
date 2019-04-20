#include <powannier.h>

using namespace POWannier;

double honeycombPotential(Position r, double Vxbar, double Vx, double Vy,
    double Vz, double kl, double theta, double phi, double alpha) {
  double x = r[0];
  double y = r[1];
  double z = r[2];
  return -Vxbar * std::pow(std::cos(kl * x + theta/2.0), 2)
          -Vx * std::pow(std::cos(kl * x), 2)
          -Vy * std::pow(std::cos(2 * kl * y), 2)
          -Vz * std::pow(std::cos(kl * z), 2)
          -2 * alpha * std::sqrt(Vx * Vy) * std::cos(kl * x) *
            std::cos(kl * y) * std::cos(phi);
}
    

int main() {
  double kl = pi;
  std::vector<Vector> basis({
    {pi/kl, -pi/kl},
    {pi/kl,  pi/kl}
  });

  auto&& func = [=] (Position r) {
    return honeycombPotential(r, 14, 0.79, 6.45, 30, 1, pi, 0.9, 0);
  };

  auto V = POWannier::Potential::create()
    .setBasis(basis)
    .evaluateFromFunction(func)
    .setCutoff(2)
    .setPrecision(1e-13)
    .complete();

  auto bs = std::make_shared<POWannier::BlochSystem>(
    POWannier::BlochSystem::createWithPotential(V)
      .setLaserWavelength(kl)
      .setFourierCutoff(6)
      .setSystemSize(10)
      .setPotentialDepth(1)
      .complete()
  );

  Wannier wan1, wan2;

  {
    POWannier::RSystem rs(bs, {0, 1});

    wan1 = rs.getWannier({0, 0}, {0, 0});
    wan2 = rs.getWannier({0, 0}, {1, 0});
  }

  POWannier::interaction(wan1, wan1, wan1, wan1);

  Position beg = {-1, -1};
  Position end = {1, 1};
  wan1.draw("wan1.dat", 100, beg, end);
  wan2.draw("wan2.dat", 100, beg, end);

  return 0;
}
