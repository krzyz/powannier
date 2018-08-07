#ifndef POWANNIER_RSYSTEM_H
#define POWANNIER_RSYSTEM_H

#include <armadillo>
#include <memory>
#include "blochsystem.h"
#include "wannier.h"

namespace POWannier {
  class RSystem {
    public:
      const int N;
      const int dim;
      const int cutoff;

      RSystem(std::shared_ptr<BlochSystem> bs);
      void generateRs();
      Wannier getWannier(NPoint n);
    private:
      std::vector<arma::vec> _riEigenvalues;
      std::vector<arma::cx_mat> _riEigenvectors;
      std::shared_ptr<BlochSystem> _bs;

      void generateR1();
      void generateSubR(arma::cx_mat eigenvalues,
                        int innerDim);

      void loopThroughMN(int innerDim, auto&& func);
  };
}

#endif
