#include "potential.h"
#include "potentialbuilder.h"

namespace POWannier {
  PotentialBuilder Potential::createWithBasis(std::vector<POWannier::Vector> basis) {
    return PotentialBuilder(basis);
  }

  double Potential::operator()(Position r) const {
    arma::cx_double result = 0;
    for (auto coefficient : _coefficients) {
      NPoint n = coefficient.first;
      arma::cx_double value = coefficient.second;
      ReciprocalVector g = _reciprocalBasis * n;
      result += value * std::exp(arma::cx_double(0, arma::dot(g, r)));
    }

    return std::real(result);
  }

  const int Potential::dim() const {
    return _dim;
  }

  const LatticeBasis& Potential::latticeBasis() const {
    return _latticeBasis;
  }

  const ReciprocalBasis& Potential::reciprocalBasis() const {
    return _reciprocalBasis;
  }

  const FourierCoefficients& Potential::fourierCoefficients() const {
    return _coefficients;
  }

  double Potential::elementaryCellVolume() const {
    return abs(arma::det(latticeBasis()));
  }

}