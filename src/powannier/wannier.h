#ifndef POWANNIER_WANNIER_H
#define POWANNIER_WANNIER_H

#include <memory>
#include "aliases.h"

namespace POWannier {
  class Wannier {
    public:
      friend double tunneling(const POWannier& wan1, const POWannier& wan2);

      virtual double operator()(Position r) const = 0;

      virtual arma::cx_double blochCoefficient(int mi, int bandi, int ni) const = 0;

      virtual void draw(std::string filename, int density, Position beg, Position end);

      virtual void draw(std::string filename, int density, DrawRange drawRange);

      virtual const std::vector<int>& bands() const = 0;

      virtual const Wannier& get1dWannier(int n) const = 0;

      virtual arma::vec energies(NPoint m) const = 0;

      virtual std::vector<POWannier::NPoint> ms() const = 0;

      int cutoff() const override;

      virtual int N() const = 0;

      virtual int dim() const = 0;
      }
  };
}
#endif


auto fun;
auto bs;

auto wan1d = SimpleWannier::createShared(bs, 5);
auto wan = WannierProduct::create(wan1d, wan1d, wan1d);

double u = interaction(wan, wan, wan, wan, 0, 0, 0, 0);