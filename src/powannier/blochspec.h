#ifndef POWANNIER_BLOCHSPEC_H
#define POWANNIER_BLOCHSPEC_H
#include <tuple>
#include "builders/potentialbuilders.h"

namespace POWannier {
  /// @cond HIDDEN
  class BlochSpec{
    public:
      const double kl;
      const int cutoff;
      const int dim;
      const double s;
      const std::vector<POWannier::NPoint> ns;
      const LatticeBasis latticeBasis;
      const ReciprocalBasis reciprocalBasis;

      BlochSpec(const Potential& V, 
          double kl, int cutoff, double s = 1);

      std::tuple<arma::vec, arma::cx_mat>
          generate(const Potential& V, ReciprocalVector k) const;

      std::complex<double> blochC(const arma::cx_mat& eigenvectors,
          NPoint n, int band) const;

      std::complex<double> bloch(const arma::cx_mat& eigenvectors,
          ReciprocalVector k, Position r, int band) const;
  };
  /// @endcond
}

#endif