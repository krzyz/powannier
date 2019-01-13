#ifndef POWANNIER_RSYSTEM_H
#define POWANNIER_RSYSTEM_H

#include <armadillo>
#include <memory>
#include "blochsystem.h"
#include "wannier.h"
#include "wannierpositions.h"
#include "../util/cache.h"

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

      Wannier getWannier(NPoint cellLocation);
      Wannier getWannier(NPoint cellLocation, NPoint wanLocation);

      const Eigensystem& getEigensystem(int n);
      const Eigensystem& getSubEigensystem(int n, NPoint cellLocation, NPoint wanLocation);

      const arma::vec& getEigenvalues(int n) {
        return std::get<0>(getEigensystem(n));
      }

      const arma::cx_mat& getEigenvectors(int n) {
        return std::get<1>(getEigensystem(n));
      }

      const arma::vec& getSubEigenvalues(int n, NPoint cellLocation, NPoint wanLocation) {
        return std::get<0>(getSubEigensystem(n, cellLocation, wanLocation));
      }

      const arma::cx_mat& getSubEigenvectors(int n, NPoint cellLocation, NPoint wanLocation) {
        return std::get<1>(getSubEigensystem(n, cellLocation, wanLocation));
      }

    private:
      // here std::vector<int> is used instead of NPoint (i.e. arma::ivec),
      // because armadillo vectors do not play well with being map keys as 
      // a part of tuple (operator< has to return bool but when comparing
      // arma vectors it rather returns a vector of bools)
      using Ipp = std::tuple<int,std::vector<int>,std::vector<int>>;

      std::shared_ptr<BlochSystem> _bs;
      WannierPositions _wannierPositions;

      // caches for storing already computed eigensystems
      Cache<int,Eigensystem> _eigensystemCache;
      Cache<Ipp,Eigensystem> _subEigensystemCache;
      Cache<Ipp,arma::cx_mat> _subspaceCache;

      arma::uvec transformToBm() const;
      // used to calculate permutations transforming bases used in calculating
      // position operator of dimension "inner" to the main basis
      arma::uvec transformFromInner(int inner) const;
      // same as above, but including the bands in the aggregate basis index
      arma::uvec transformFromInnerToBm(int inner) const;

      void sortEigensystem(arma::vec& eigenvalues, arma::cx_mat& eigenvectors) const;

      Eigensystem calculateEigensystem(int inner) const;
      Eigensystem calculateSubEigensystem(int inner, const arma::cx_mat& subspace);
      Eigensystem calculateSubEigensystem(Ipp ipp);
      arma::cx_mat calculateSubspace(Ipp ipp);
      std::tuple<int,int> getSubspaceRange(int n, NPoint cellLocation, NPoint wanLocation, bool last = false);
      const arma::cx_mat& getSubspace(int n, NPoint cellLocation, NPoint wanLocation);
      int wannierPosition(int n, NPoint cellLocation, NPoint wanLocation);

      // calculate part of the block in the position operator of dimension "inner"
      //   block corresponds to m indices (m_{i} for i = inner)
      //   part corresponds to bands bi (row) and bj (column)
      // m are integers that reproduce wavewectors from first brillouin zone:
      //   k = (2 * \pi * (m - f(N)) / N, where f(n) = 
      arma::cx_mat rInnerMatrix(int inner, NPoint m, int bi, int bj) const;

      // summation over all fourier coefficents for indices m (row) and mp (column),
      // and bands bands[bi] (row) and bands[bj] (column)
      arma::cx_double sumOverNs(int inner, NPoint m, NPoint mp, int bi, int bj) const;




      arma::imat basisMs(int inner) const;
  };
}

#endif
