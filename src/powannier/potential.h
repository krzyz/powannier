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

  template <class IntegrationProvider>
  class PotentialT {
    public:
      const int dim;

      template<class Function>
      PotentialT(
          std::vector<Vector> latticeBasis,
          Function&& potentialFunction,
          int cutoff, double prec = 1e-15);

      double operator()(Position r);
      double relative_error(double prec = 1e-14);

    private:
      double _prec;
      LatticeBasis _latticeBasis;
      ReciprocalBasis _reciprocalBasis;
      FourierCoefficients _potentialCoefficients;
      std::function<double(Position)> _potentialFunction;
      std::map<std::string, double> _cache;

      template<class Function>
      FourierCoefficients calculateFourierCoefficients(
          Function&& potentialFunction, int cutoff);
  };

  template <class IntegrationProvider>
  template <class Function>
  PotentialT<IntegrationProvider>::PotentialT(
      std::vector<Vector> latticeBasis,
      Function&& potentialFunction,
      int cutoff, double prec) :
    dim(latticeBasis.size()),
    _prec(prec),
    _potentialFunction(potentialFunction) {
    if (std::is_base_of<IntegrationInterface, IntegrationProvider>::value == false) {
      throw std::logic_error("IntegrationProvider template argument must implement POWannier::IntegrationInterface!");
    }

    std::vector<double> basisFlattened;
    for (auto rvec : latticeBasis) {
      for (auto val : rvec) {
        basisFlattened.push_back(val);
      }
    }

    _latticeBasis = LatticeBasis(&basisFlattened.front(), dim, dim).t();
    _reciprocalBasis = 2 * pi * arma::inv(_latticeBasis);

    _potentialCoefficients = calculateFourierCoefficients(potentialFunction, cutoff);
  }

  template <class IntegrationProvider>
  template <class Function>
  FourierCoefficients PotentialT<IntegrationProvider>::calculateFourierCoefficients(
      Function&& potentialFunction, int cutoff) {
    FourierCoefficients coefficients;
    Position xmin(dim, arma::fill::zeros);
    Position xmax(dim, arma::fill::ones);

    for (NPoint n : nspace(cutoff, dim)) {
      double i = 0;
      arma::cx_double coefficient = chop(IntegrationProvider::integrate(
          [&, n] (Position alpha) {
            Position r = alpha * _latticeBasis;
            return potentialFunction(r) * 1.0 *
                std::exp(arma::cx_double(0, -2 * pi * arma::dot(alpha, n)));
          }, xmin, xmax, _prec, _prec));

      if (coefficient != 0.0) {
        coefficients.push_back(std::make_pair(n, coefficient));
      }
    }

    return coefficients;
  }


  template <class IntegrationProvider>
  double PotentialT<IntegrationProvider>::operator()(Position r) {
    arma::cx_double result = 0;
    for (auto coefficient : _potentialCoefficients) {
      NPoint n = coefficient.first;
      arma::cx_double value = coefficient.second;
      ReciprocalVector g = _reciprocalBasis * n;
      result += value * std::exp(arma::cx_double(0, arma::dot(g, r)));
    }
    return std::real(result);
  }


  template <class IntegrationProvider>
  double PotentialT<IntegrationProvider>::relative_error(double prec) {
    if (_cache.find("relative_error") != _cache.end()) {
      return _cache["relative_error"];
    } else {
      Position xmin(dim, arma::fill::zeros);
      Position xmax(dim, arma::fill::ones);

      double difference_squared = IntegrationProvider::integrate([&] (Position alpha) {
          Position r = alpha * _latticeBasis;
          return std::pow(_potentialFunction(r) - (*this)(r), 2);
          }, xmin, xmax, prec, prec);

      double i = 0;
      double norm = IntegrationProvider::integrate([&] (Position alpha) {
          ++i;
          Position r = alpha * _latticeBasis;
          return std::pow(_potentialFunction(r), 2);
          }, xmin, xmax, prec, prec);

      double relative_error = difference_squared / norm;

      _cache["relative_error"] = relative_error;

      return relative_error;
    }
  }


  using Potential = PotentialT<CubatureIntegration>;
}

#endif
