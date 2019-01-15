#include "potentialfunctionbuilder.h"

namespace POWannier {
  PotentialFunctionBuilder& PotentialFunctionBuilder::setCutoff(int cutoff) {
    _cutoff = cutoff;
    return *this;
  }

  PotentialFunctionBuilder& PotentialFunctionBuilder::setPrecision(double precision) {
    _precision = precision;
    return *this;
  }

  PotentialFunctionBuilder& PotentialFunctionBuilder::setRelativeErrorPrecision(double precision) {
    _relErrorPrecision = precision;
    return *this;
  }


  PotentialFunctionBuilder& PotentialFunctionBuilder::requireRelativeErrorBelow(double maxError) {
    _maxError = maxError;
    return *this;
  }
}