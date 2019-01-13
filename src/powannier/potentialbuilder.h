#ifndef POWANNIER_POTENTIALBUILDER_H
#define POWANNIER_POTENTIALBUILDER_H

#include <vector>
#include <queue>

#include "potential.h"

namespace POWannier {
  class Potential;
  class PotentialFourierBuilder;
  class PotentialFunctionBuilder;

  class PotentialBuilder {
    public:
      PotentialBuilder();

      PotentialBuilder& setBasis(std::vector<POWannier::Vector> basis);
      PotentialFourierBuilder addFourierCoefficients();

      template <class Function>
      PotentialFunctionBuilder evaluateFromFunction(Function&& func);
      
      Potential&& complete();

    protected:
      Potential& _potential;
      explicit PotentialBuilder(Potential& potential) :
        _potential(potential) {}

    private:
      enum class State {
        BasisNotDeclared,
        BasisDeclared,
      };

      void ensureBasisDefined();

      Potential _p;
      State _state;
  };

  class PotentialFourierBuilderAddition;

  class PotentialFourierBuilder : private PotentialBuilder {
    public:
      explicit PotentialFourierBuilder(Potential& potential) :
        PotentialBuilder(potential) {
      }

      PotentialFourierBuilderAddition provideMultiple(std::vector<NPoint> indices);
      PotentialFourierBuilderAddition provideOne(NPoint index);

      using PotentialBuilder::complete;
  };

  class PotentialFourierBuilderAddition : private PotentialBuilder {
    public:
      explicit PotentialFourierBuilderAddition(
        Potential& potential, std::queue<NPoint> indices) :
        PotentialBuilder(potential), _indices(indices) {}
      
      PotentialFourierBuilder values(std::vector<Complex> values);
      PotentialFourierBuilder value(Complex value);

    private:
      void addToPotentialCoefficients(Complex value);
      std::queue<NPoint> _indices;
  };

  class PotentialFunctionBuilder : private PotentialBuilder {
    public:
      template <class Function>
      explicit PotentialFunctionBuilder(Potential& potential, Function&& func) :
        PotentialBuilder(potential),
        _potentialFunction(func) {}

      PotentialFunctionBuilder& setCutoff(int cutoff);
      PotentialFunctionBuilder& setPrecision(double precision);
      PotentialFunctionBuilder& setRelativeErrorPrecision(double precision);
      PotentialFunctionBuilder& requireRelativeErrorBelow(double maxError);

      template <class IntegrationProvider = CubatureIntegration>
      double relativeError();
     
      template <class IntegrationProvider = CubatureIntegration>
      Potential&& complete();


    private:
      std::function<double(Position)> _potentialFunction;
      int _cutoff = 1;
      double _precision = 1e-15;
      double _relErrorPrecision = 1e-14;
      double _maxError = -1;
  };

  template <class Function>
  PotentialFunctionBuilder PotentialBuilder::evaluateFromFunction(Function&& func) {
    ensureBasisDefined();
    return PotentialFunctionBuilder(_potential, func);
  }

  template <class IntegrationProvider>
  double PotentialFunctionBuilder::relativeError() {
    Position xmin(_potential.dim(), arma::fill::zeros);
    Position xmax(_potential.dim(), arma::fill::ones);

    double difference_squared = IntegrationProvider::integrate([&] (Position alpha) {
        Position r = alpha * _potential.latticeBasis();
        return std::pow(_potentialFunction(r) - _potential(r), 2);
        }, xmin, xmax, _relErrorPrecision, _relErrorPrecision);

    double i = 0;
    double norm = IntegrationProvider::integrate([&] (Position alpha) {
        ++i;
        Position r = alpha * _potential.latticeBasis();
        return std::pow(_potentialFunction(r), 2);
        }, xmin, xmax, _relErrorPrecision, _relErrorPrecision);

    double relative_error = difference_squared / norm;

    return relative_error;
  }

  template <class IntegrationProvider>
  Potential&& PotentialFunctionBuilder::complete() {
    _potential._coefficients.clear();
    Position xmin(_potential.dim());
    Position xmax(_potential.dim());
    xmin.fill(-0.5);
    xmax.fill(0.5);

    for (NPoint n : nspace(_cutoff, _potential.dim())) {
      double i = 0;
      arma::cx_double coefficient = chop(IntegrationProvider::integrate(
          [&, n] (Position alpha) {
            Position r = alpha * _potential.latticeBasis();
            return _potentialFunction(r) * 1.0 *
                std::exp(arma::cx_double(0, -2 * pi * arma::dot(alpha, n)));
          }, xmin, xmax, _precision, _precision));

      if (coefficient != 0.0) {
        _potential._coefficients.push_back(std::make_pair(n, coefficient));
      }
    }

    if (_maxError >= 0) {
      if (relativeError<IntegrationProvider>() > _maxError) {
        throw std::runtime_error("Relative error of approximated potential exceeded the requirements.");
      }
    }

    return std::move(PotentialBuilder::complete());
  }
}

#endif
