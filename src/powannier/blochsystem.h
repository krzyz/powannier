#ifndef POWANNIER_BLOCHSYSTEM_H
#define POWANNIER_BLOCHSYSTEM_H
#include <memory>
#include <mutex>
#include "potential.h"

namespace POWannier { 
  class BlochSystem {
    public:
      const double kl;
      const int cutoff;
      const int N;
      const int dim;
      const double s;
      std::vector<POWannier::NPoint> ms;
      std::vector<POWannier::NPoint> ns;
      std::shared_ptr<const Potential> V;

      BlochSystem(std::shared_ptr<const Potential> V, double kl, int cutoff, int N, double s = 1);
      void generate(NPoint m);
      void generateAll();
      const arma::vec& energies(NPoint m);
      const arma::cx_mat& eigenvectors(NPoint m);
      std::complex<double> bloch(NPoint m, Position r, int band = 0);
      std::complex<double> blochC(NPoint m, NPoint n, int band = 0);
      ReciprocalVector kFromM(NPoint m);

    private:
      std::vector<arma::vec> _energies;
      std::vector<arma::cx_mat> _eigenvectors;
      ReciprocalBasis _reciprocalBasis;

      arma::cx_mat generateMatrix(const Potential& V, double s = 1);
  };
}

#endif
