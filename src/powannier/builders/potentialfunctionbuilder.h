#ifndef _POWANNIER_POTENTIALFUNCTIONBUILDER_H
#define _POWANNIER_POTENTIALFUNCTIONBUILDER_H

#include "potentialbuilder.h"

namespace POWannier {
  /**
   * @brief Builder of the Potential class by integrating the provided
   * function
   */
  class PotentialFunctionBuilder : private PotentialBuilder {
    public:
      /** @brief Exception thrown when requirements for relative error of 
       * the potential approximation are not met.
       */
      class RelativeErrorNotMetError : public std::runtime_error {
        public:
          /// Constructor.
          RelativeErrorNotMetError() : 
            std::runtime_error("Relative error requirements not met!") {}
      };
      /**
       * @brief Basic setter.
       * The default value is 1;
       * @param cutoff
       *   Cutoff in plane wave expansion of the potential (see Bloch).
       * @returns
       *   Reference to the current builder.
       */
      PotentialFunctionBuilder& setCutoff(int cutoff);

      /**
       * @brief Basic setter.
       * The default value is \f$ 10^{-15} \f$.
       * @param precision
       *   Precision used in the integration.
       * @returns
       *   Reference to the current builder.
       */
      PotentialFunctionBuilder& setPrecision(double precision);

      /**
       * @brief Basic setter.
       * The default value is \f$ 10^{-14} \f$.
       * @param precision
       *   Precision used in calculating the relative error of
       *   the final potential approximation. 
       * @returns
       *   Reference to the current builder.
       */
      PotentialFunctionBuilder& setRelativeErrorPrecision(double precision);

      /**
       * @brief Basic setter.
       * By default, no requirements are set.
       * @param maxError
       *   A value over which 
       * @returns
       *   Reference to the current builder.
       */
      PotentialFunctionBuilder& requireRelativeErrorBelow(double maxError);

      /**
       * @brief Finish building the Potential and return the result.
       * 
       * After the call to this method, the Fourier coefficients for the 
       * potential are calculated using @p IntegrationProvider. The cutoff
       * and integration precision can be set by calling setter methods 
       * before calling this function. \n 
       * 
       * If function requreRelativeErrorBelow() was used, then after
       * calculating Fourier coefficients a relative error is calculated
       * (using the @p IntegrationProvider class) and if the requirements
       * are not met, the RelativeErrorNotMetError is thrown.

       * @tparam IntegrationProvider
       *   Class used for integration
       * 
       * @returns
       *   Finished object of Potential class.
       * 
       * @pre
       *   @p IntegrationProvider class has to implement 
       *   IntegrationProvider::integrate() method
       * 
       * @throws
       *   RelativeErrorNotMetError
       */
      template <class IntegrationProvider = CubatureIntegration>
      Potential&& complete();

    private:
      friend class PotentialBuilder;

      template <class Function>
      explicit PotentialFunctionBuilder(Potential& potential, Function&& func) :
        PotentialBuilder(potential),
        _potentialFunction(func) {}

      template <class IntegrationProvider = CubatureIntegration>
      double relativeError();

      std::function<double(Position)> _potentialFunction;
      int _cutoff = 1;
      double _precision = 1e-15;
      double _relErrorPrecision = 1e-14;
      double _maxError = -1;
  };

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

    // calculate Fourier series coefficients (V_n) for a given potential,
    // that is V(r) = \sum_n V_n e^{i (\sum_i n_i \vec{b_i}) \vec{r}},
    // where n is a set of indices (summation is performed for all 
    // integer indices in the range [-cutoff, cutoff]), 
    // and \vec{b_i} are reciprocal lattice basis vectors
    for (NPoint n : nspace(_cutoff, _potential.dim())) {
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

  template <class Function>
  PotentialFunctionBuilder PotentialBuilder::evaluateFromFunction(Function&& func) {
    ensureBasisDefined();
    return PotentialFunctionBuilder(_potential, func);
  }
}

#endif