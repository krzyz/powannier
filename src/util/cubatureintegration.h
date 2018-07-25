#ifndef POWANNIER_CUBATUREINTEGRATION_H
#define POWANNIER_CUBATUREINTEGRATION_H

#include <armadillo>
#include <complex>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "integrationinterface.h"
#include "cubature.h"

namespace POWannier {
  /// @cond HIDDEN
  template<class T> struct is_complex : std::false_type {};
  template<class T> struct is_complex<std::complex<T>> : std::true_type {};
  template<class T> struct is_cont : std::false_type {};
  template<class T, class Alloc> struct is_cont<std::vector<T, Alloc>> : std::true_type {};
  /// @endcond

  /**
   * \brief An adapter implementing IntegrationInterface using Cubature C library (https://github.com/stevengj/cubature).
   */

  class CubatureIntegration : IntegrationInterface {
    public:
      template <class Function, class Vector, typename std::enable_if_t<is_complex<typename std::result_of<Function(Vector)>::type>{}>* = nullptr>
      static typename std::result_of<Function(Vector)>::type integrate(Function&& function, Vector xmin, Vector xmax) {
        using namespace std::literals;

        double resultReal, resultImag, err;

        auto functionReal = [&] (Vector x) -> double {
            return std::real(function(x));
        };

        auto functionImag = [&] (Vector x) -> double {
            return std::imag(function(x));
        };

        int dim = checkLimits(xmin, xmax);

        pcubature(1, integrateForward<decltype(functionReal), Vector>, &functionReal, dim, toArray(xmin), toArray(xmax), 0, 1e-15, 1e-15, ERROR_INDIVIDUAL, &resultReal, &err);
        if (abs(resultReal) < 1e-14) {
          resultReal = 0;
        }
        pcubature(1, integrateForward<decltype(functionImag), Vector>, &functionImag, dim, toArray(xmin), toArray(xmax), 0, 1e-15, 1e-15, ERROR_INDIVIDUAL, &resultImag, &err);
        if (abs(resultImag) < 1e-14) {
          resultReal = 0;
        }

        return resultReal + resultImag * 1i;
      }

      template <class Function, class Vector, typename std::enable_if_t<!is_complex<typename std::result_of<Function(Vector)>::type>{}>* = nullptr>
      static typename std::result_of<Function(Vector)>::type integrate(Function&& function, Vector xmin, Vector xmax) {
        int dim = checkLimits(xmin, xmax);

        double result, err;
        pcubature(1, integrateForward<Function, Vector>, &function, dim, toArray(xmin), toArray(xmax), 0, 1e-15, 1e-15, ERROR_INDIVIDUAL, &result, &err);

        return result;
      }


    private:
      template <class Vector, typename std::enable_if_t<is_cont<Vector>{}>* = nullptr>
      static double* toArray(Vector x) {
          return x.data();
      }

      template <class Vector, typename std::enable_if_t<!is_cont<Vector>{}>* = nullptr>
      static double* toArray(Vector x) {
          return x.memptr();
      }

      template <class Vector, typename std::enable_if_t<is_cont<Vector>{}>* = nullptr>
      static int checkLimits(Vector xmin, Vector xmax) {
        int dim = xmin.size();
        if (xmax.size() != xmin.size()) {
          throw std::invalid_argument("Dimensions of integral boundaries must be equal!");
        }

        return dim;
      }

      template <class Vector, typename std::enable_if_t<!is_cont<Vector>{}>* = nullptr>
      static int checkLimits(Vector xmin, Vector xmax) {
        int dim = xmin.n_elem;
        if (arma::size(xmax) != arma::size(xmax)) {
          throw std::invalid_argument("Dimensions of integral boundaries must be equal!");
        }

        return dim;
      }


      template <class Function, class Vector>
      static int integrateForward(unsigned ndim, const double* xArray, void* fdata, unsigned fdim, double* fval);
  };

  template <class Function, class Vector>
  int CubatureIntegration::integrateForward(unsigned ndim, const double* xArray, void* fdata, unsigned fdim, double* fval) {
    using Base_Vector = typename std::remove_cv<typename std::remove_reference<Vector>::type>::type;

    Base_Vector x(ndim);

    for (int i = 0; i < ndim; ++i) {
      x[i] = xArray[i];
    }

    *fval = (*reinterpret_cast<Function*&&>(fdata))(x);
    return 0;
  }
}

#endif
