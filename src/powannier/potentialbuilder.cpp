#include "potentialbuilder.h"

namespace POWannier {
  PotentialBuilder::PotentialBuilder() : 
    _potential(_p), 
    _state(State::BasisNotDeclared) {
    }

  PotentialBuilder& PotentialBuilder::setBasis(std::vector<POWannier::Vector> basis) {
    _potential._dim = basis.size();
    std::vector<double> basisFlattened;
    for (auto rvec : basis) {
      for (auto val : rvec) {
        basisFlattened.push_back(val);
      }
    }

    _potential._latticeBasis = LatticeBasis(
      &basisFlattened.front(), _potential.dim(), _potential.dim()).t();
    _potential._reciprocalBasis = 2 * pi * arma::inv(_potential._latticeBasis);

    _state = State::BasisDeclared;

    return *this;
  }

  PotentialFourierBuilder PotentialBuilder::addFourierCoefficients() {
    ensureBasisDefined();
    return PotentialFourierBuilder{_potential};
  }

  Potential&& PotentialBuilder::complete() {
    return std::move(_potential);
  }

  void PotentialBuilder::ensureBasisDefined() {
    if (_state != State::BasisDeclared) {
      throw std::runtime_error("Attempted to take action that requires a basis without supplying one");
    }
  }

  PotentialFourierBuilderAddition PotentialFourierBuilder::provideMultiple(std::vector<NPoint> indices) {
    std::queue<NPoint> indicesQ;

    for (auto index : indices) {
      indicesQ.push(index);
    }

    return PotentialFourierBuilderAddition{_potential, indicesQ};
  }

  PotentialFourierBuilderAddition PotentialFourierBuilder::provideOne(NPoint index) {
    std::queue<NPoint> indices;
    indices.push(index);
    return PotentialFourierBuilderAddition{_potential, indices};
  }

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