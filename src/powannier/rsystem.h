#ifndef POWANNIER_RSYSTEM_H
#define POWANNIER_RSYSTEM_H

#include <armadillo>
#include <memory>
#include "blochsystem.h"
#include "wannier.h"
#include "wannierpositions.h"

namespace POWannier {
  class RSystem {
    public:
      const int N;
      const int dim;
      const int mdim;
      const int bmdim;
      const std::vector<int> bands;
      const int cutoff;

      RSystem(std::shared_ptr<BlochSystem> bs,
        std::vector<int> bands);

      RSystem(std::shared_ptr<BlochSystem> bs,
        std::vector<int> bands,
        WannierPositions positions);
      Wannier getWannier(NPoint n);
      Wannier getWannier(NPoint n, NPoint positions);
      const arma::cx_mat& r1Eigenvectors();
      const arma::vec& r1Eigenvalues();

    private:
      EigenSystem _r1EigenSystem;
      std::shared_ptr<BlochSystem> _bs;
      WannierPositions _wannierPositions;
      arma::uvec transformToBm() const;
      arma::uvec transformFromInner(int inner) const;
      arma::uvec transformFromInnerToBm(int inner) const;

      void sortEigensystem(arma::vec& eigenvalues, arma::cx_mat& eigenvectors) const;
      EigenSystem getREigensystem(int inner) const;
      EigenSystem getSubREigensystem(int inner, const arma::cx_mat& subspace) const;
      arma::cx_mat rInnerMatrix(int inner, NPoint m, int bi, int bj) const;
      arma::cx_double sumOverNs(int inner, NPoint m, NPoint np, int bi, int bj) const;
  };
}

#endif
