#include <armadillo>
#include <exception>
#include <utility>
#include <vector>

#include "aliases.h"
#include "helpers.h"

template <class IntegrationProvider>
namespace POWannier {

  template <class IntegrationProvider>
  class PotentialT {
    public:
      template<class Function>
      PotentialT(
          std::vector<Vector> latticeBasis,
          Function&& potentialFunction,
          int cutoff);

      double operator()(Position r);

    private:
      const int dim;
      LatticeBasis _latticeBasis;
      ReciprocalBasis _reciprocalBasis;
      FourierCoefficients _potentialCoefficients;
      std::function<double(Position)> _potentialFunction;

      template<class Function>
      calculateFourierCoefficients(
          Function&& potentialFunction, int cutoff);
  };

  template <class IntegrationProvider>
  template <class Function>
  PotentialT::PotentialT(
      std::vector<Vector> latticeBasis,
      Function&& potentialFunction,
      int cutoff) :
    dim(latticeBasis.size()),
    _potentialFunction(potentialFunction) {
    if (std::is_base_of<POWannier::IntegrationInterface, IntegrationProvider>::value == false) {
      throw std::logic_error("IntegrationProvider template argument must implement POWannier::IntegrationInterface!");
    }
    for (auto rvec : latticeBasis) {
      for (auto val : rvec) {
        basisFlattened.push_back(val);
      }
    }

    _latticeBasis = LatticeBasis(&basisFlattened.front(), dim, dim);
    _reciprocalBasis = 2 * pi * arma::inv(_latticeBasis);

    _potentialCoefficients = calculateFourierCoefficients(potentialFunction, cutoff);
  }

  template <class IntegrationProvider>
  template <class Function>
  FourierCoefficients PotentialT::calculateFourierCoefficients(
      Function&& potentialFunction, int cutoff) {
    using namespace std::string_literals;
    FourierCoefficients coefficients;
    Position xmin(dim, arma::fill::zeros);
    Position xmax(dim, arma::fill::ones);

    for (NPoint n : nspace(cutoff, dim)) {
      arma::cx_double coefficient = chop(IntegrationProvider::integrate(
          [this, &, n] (Position alpha) {
            Position r = alpha * _latticeBasis;
            return potentialFunction(r) * 
                exp(0, 1i * 2 * pi * arma::dot(alpha, n));
          }, xmin, xmax));

      if (coefficient != 0.0) {
        coefficients.push_back(std::make_pair(n, coefficient));
      }
    }

    return coefficients;
  }


  template <class IntegrationProvider>
  double PotentialT::operator()(Position r) {
    using namespace std::string_literals;
    double result = 0;
    for (auto coefficient : coefficients) {
      NPoint n = coefficient.first;
      arma::cx_double value = coefficient.second;
      ReciprocalVector g = _reciprocalBasis * n;
      result += coefficient * exp(-1i * arma::dot(g, r));
    }
    return result;
  }
}
