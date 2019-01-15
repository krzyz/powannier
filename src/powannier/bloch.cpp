#include "bloch.h"

namespace POWannier {
  Bloch::Bloch(ReciprocalVector k, const Potential& V, double laserWavelength,
      int fourierCutoff, double potentialDepth) :
    _blochSpec(V, laserWavelength, fourierCutoff, potentialDepth),
    _k(k) {
    std::tie(_energies, _eigenvectors) = _blochSpec.generate(V, k);
  }

  std::complex<double> Bloch::operator()(Position r, int band) {
    return _blochSpec.bloch(_eigenvectors, _k, r, band);
  }
}