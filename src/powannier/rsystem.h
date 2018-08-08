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
      Wannier getWannier(NPoint n);
      const arma::cx_mat& r1Eigenvectors();
      const arma::vec& r1Eigenvalues();

    private:
      EigenSystem _r1EigenSystem;
      std::shared_ptr<BlochSystem> _bs;

      void sortEigensystem(arma::vec& eigenvalues, arma::cx_mat& eigenvectors);
      arma::uvec rearrangedIndicesR1();
      EigenSystem generateR1();
      EigenSystem generateSubR(arma::cx_mat eigenvalues,
                        int innerDim);

      void loopThroughMN(int innerDim, auto&& func);
      //void test();
  };
}

#endif
