#include "bloch.h"

namespace POWannier {
  Bloch::Bloch(ReciprocalVector k, const Potential& V, double kl, int cutoff, double s) :
    _blochSpec(V, kl, cutoff, s),
    _k(k) {
    std::tie(_energies, _eigenvectors) = _blochSpec.generate(V, k);
  }

  std::complex<double> Bloch::bloch(Position r, int band) {
    return _blochSpec.bloch(_eigenvectors, _k, r, band);
  }
}