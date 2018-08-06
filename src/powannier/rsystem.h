#ifndef POWANNIER_RSYSTEM_H
#define POWANNIER_RSYSTEM_H

#include <armadillo>
#include <memory>
#include "blochsystem.h"

namespace POWannier {
  class RSystem {
    public:
      const int N;
      const int dim;
      const int cutoff;

      RSystem(std::shared_ptr<BlochSystem> bs);
      void generateR();
    private:
      std::vector<arma::sp_cx_mat> _ris;
      std::shared_ptr<BlochSystem> _bs;

      void loopThroughMN(int innerDim, auto&& func);
  };
}

#endif
