#include "blochsystem.h"
#include "helpers.h"


namespace POWannier {
  BlochSystem::BlochSystem(const Potential& V, double kl, int cutoff, int N, double s) :
    N(N),
    ms(mspace(N, V.dim())),
    _elementaryCellVolume(V.elementaryCellVolume()),
    _latticeBasis(V.latticeBasis()),
    _blochSpec(V, kl, cutoff, s) {
      generateAll(V);
    }

  BlochSystemBuilder BlochSystem::createWithPotential(const Potential& potential) {
    return BlochSystemBuilder(potential);
  }

  void BlochSystem::generateAll(const Potential& potential) {
    _energies = std::vector<arma::vec>(ms.size());
    _eigenvectors = std::vector<arma::cx_mat>(ms.size());

    for (size_t mi = 0; mi < ms.size(); ++mi) {
      const auto& m = ms[mi];
      ReciprocalVector k = kFromM(m);

      std::tie(_energies[mi], _eigenvectors[mi]) = _blochSpec.generate(potential, k);
    }
  }

  const arma::vec& BlochSystem::energies(NPoint m) const {
    return _energies[mIndex(m, N)];
  }
  
  const arma::cx_mat& BlochSystem::eigenvectors(NPoint m) const {
    return _eigenvectors[mIndex(m, N)];
  }

  std::complex<double> BlochSystem::blochC(NPoint m, NPoint n, int band) const {
    return _blochSpec.blochC(eigenvectors(m), n, band);
  }

  std::complex<double> BlochSystem::bloch(NPoint m, Position r, int band) const {
    ReciprocalVector k = kFromM(m);
    return _blochSpec.bloch(eigenvectors(m), k, r, band);
  }

  ReciprocalVector BlochSystem::kFromM(NPoint m) const {
    return _blochSpec.reciprocalBasis * m / N;
  }
}
