#include "potentialbuilder.h"
#include "potentialfourierbuilder.h"
#include "potentialfunctionbuilder.h"

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
}