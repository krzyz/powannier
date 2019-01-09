#ifndef POWANNIER_POTENTIAL_H
#define POWANNIER_POTENTIAL_H

#include <armadillo>
#include <exception>
#include <map>
#include <utility>
#include <vector>

#include "../util/integrationinterface.h"
#include "../util/cubatureintegration.h"
#include "aliases.h"
#include "helpers.h"

namespace POWannier {
  class PotentialBuilder;

  class Potential {
    friend class PotentialBuilder;
    friend class PotentialFourierBuilderAddition;
    friend class PotentialFunctionBuilder;

    public:
      Potential(const Potential&) = default;
      Potential(Potential&&) = default;

      static PotentialBuilder create();
      double operator()(Position r) const;
      const int dim() const;
      const LatticeBasis& latticeBasis() const;
      const ReciprocalBasis& reciprocalBasis() const;
      const FourierCoefficients& fourierCoefficients() const;
      double elementaryCellVolume() const;

    private:
      int _dim;
      LatticeBasis _latticeBasis;
      ReciprocalBasis _reciprocalBasis;
      FourierCoefficients _coefficients;

      Potential() {}
  };
}

#endif
