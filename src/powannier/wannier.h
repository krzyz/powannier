#ifndef POWANNIER_WANNIER_H
#define POWANNIER_WANNIER_H

#include <memory>
#include "blochsystem.h"

namespace POWannier {
  enum DrawRange {
    WholeLattice,
    ElementaryCell,
  };
  class Wannier {
    public:
      const std::vector<int> bands;
      Wannier(std::shared_ptr<BlochSystem> bs, 
              const arma::cx_vec& coefficients,
              std::vector<int> bands = {0});
      double operator()(Position r);

      void drawB(std::string filename, int density, Position beg, Position end);
      void drawB(std::string filename, int density, DrawRange drawRange);
    private:
      std::shared_ptr<BlochSystem> _bs;
      arma::cx_vec _coefficients;
  };
}
#endif