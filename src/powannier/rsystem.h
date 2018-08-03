#ifndef POWANNIER_RSYSTEM_H
#define POWANNIER_RSYSTEM_H

#include <armadillo>
#include "blochsystem.h"

namespace POWannier {
  class RSystem {
    public:
    private:
      std::vector<arma::sp_cx_mat> _ris;
      void generateRHamiltonians(BlochSystem bs);
  };
}

#endif