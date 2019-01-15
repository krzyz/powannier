#include "potentialfourierbuilder.h"
#include "potentialfourierbuilderaddition.h"

namespace POWannier {
  void PotentialFourierBuilderAddition::addToPotentialCoefficients(Complex value) {
    if (value != 0.0) {
      _potential._coefficients.push_back(std::make_pair(_indices.front(), value));
    }
    _indices.pop();
  }

  PotentialFourierBuilder PotentialFourierBuilderAddition::values(std::vector<Complex> values) {
    if (_indices.size() - values.size() < 0) {
      throw std::logic_error("More Fourier coefficient values were given then indices specified.");
    }

    for (auto value : values) {
      addToPotentialCoefficients(value);
    }

    return PotentialFourierBuilder{_potential};
  }

  PotentialFourierBuilder PotentialFourierBuilderAddition::value(Complex value) {
    if (_indices.size() < 1) {
      throw std::logic_error("Attempted to give a Fourier coefficient value without specifying its index.");
    }

    addToPotentialCoefficients(value);

    return PotentialFourierBuilder{_potential};
  }
}