#ifndef POWANNIER_BLOCHSYSTEM_H
#define POWANNIER_BLOCHSYSTEM_H

#include "blochspec.h"
#include "blochsystembuilder.h"

namespace POWannier { 
  class BlochSystem {
    public:
      const int N;
      const std::vector<POWannier::NPoint> ms;

      static BlochSystemBuilder createWithPotential(const Potential& potential);

      const arma::vec& energies(NPoint m) const;
      const arma::cx_mat& eigenvectors(NPoint m) const;
      std::complex<double> bloch(NPoint m, Position r, int band = 0) const;
      std::complex<double> blochC(NPoint m, NPoint n, int band = 0) const;
      ReciprocalVector kFromM(NPoint m) const;

      double kl() const {
        return _blochSpec.kl;
      }

      int cutoff() const {
        return _blochSpec.cutoff;
      }

      int dim() const {
        return _blochSpec.dim;
      }

      double elementaryCellVolume() const {
        return _elementaryCellVolume;
      }

      const LatticeBasis& latticeBasis() const {
        return _latticeBasis;
      }

    private:
      const double _elementaryCellVolume;
      const LatticeBasis _latticeBasis;
      const BlochSpec _blochSpec;
      std::vector<arma::vec> _energies;
      std::vector<arma::cx_mat> _eigenvectors;

      void generateAll(const Potential& potential);
      BlochSystem(const Potential& potential, 
        double kl, int cutoff, int N, double s = 1);

    friend class BlochSystemBuilder;
  };
}

#endif
