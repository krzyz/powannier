#ifndef POWANNIER_BLOCH_H
#define POWANNIER_BLOCH_H
#include "blochspec.h"

namespace POWannier {
  class Bloch {
    public:
      Bloch(ReciprocalVector k, const Potential& V,
        double kl, int cutoff, double s = 1);

      ReciprocalVector k() const {
        return _k;
      }

      arma::vec energies() const {
        return _energies;
      }

      std::vector<POWannier::NPoint> ns() const {
        return _blochSpec.ns;
      };

      std::complex<double> bloch(Position r, int band);

    private:
      const BlochSpec _blochSpec;
      const ReciprocalVector _k;
      arma::vec _energies;
      arma::cx_mat _eigenvectors;
  };
}

#endif