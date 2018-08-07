#ifndef POWANNIER_WANNIER_H
#define POWANNIER_WANNIER_H

#include <memory>
#include "blochsystem.h"

namespace POWannier {
  class Wannier {
    public:
      Wannier(std::shared_ptr<BlochSystem> bs, 
              const arma::cx_vec& coefficients);
      double operator()(Position r);

    private:
      std::shared_ptr<BlochSystem> _bs;
      arma::cx_vec _coefficients;
  };
}
#endif