#ifndef POWANNIER_CUBATUREINTEGRATION_H
#define POWANNIER_CUBATUREINTEGRATION_H

#include <armadillo>
#include <complex>
#include <stdexcept>
#include <type_traits>
#include <random>
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
   * \brief Adapter implementing IntegrationInterface using Cubature C library (https://github.com/stevengj/cubature).
   */

  class CubatureIntegration : IntegrationInterface {
    public:
      /**
       * Implementation of IntegrationInterface::integrate() function for complex valued integrand.
       */
      template <class Function, class Vector, typename std::enable_if_t<is_complex<typename std::result_of<Function(Vector)>::type>{}>* = nullptr>
      static typename std::result_of<Function(Vector)>::type integrate(Function&& function, Vector xmin, Vector xmax, double prec_rel = 1e-15, double prec_abs = 1e-15, int level = 0) {
        using namespace std::literals;
        
        int dim = checkLimits(xmin, xmax);

        double resultReal, resultImag, err;

        auto functionReal = [&] (Vector x) -> double {
            return std::real(function(x));
        };

        auto functionImag = [&] (Vector x) -> double {
            return std::imag(function(x));
        };

        if (checkIfIntegrationFirstStepOk(functionReal, xmin, xmax, dim, prec_rel, prec_abs, level) == false ||
            checkIfIntegrationFirstStepOk(functionImag, xmin, xmax, dim, prec_rel, prec_abs, level) == false) {
          return splitIntegrate(std::ref(functionReal), xmin, xmax, dim, prec_rel, prec_abs, level) + 
            splitIntegrate(std::ref(functionImag), xmin, xmax, dim, prec_rel, prec_abs, level)* 1i;
        }

        pcubature_v(1, reinterpret_cast<cubatureFunctionPointerType>(integrateForward<decltype(functionReal), Vector>), &functionReal, dim, toArray(xmin), toArray(xmax), 0, prec_rel, prec_abs, ERROR_INDIVIDUAL, &resultReal, &err);

        pcubature_v(1, reinterpret_cast<cubatureFunctionPointerType>(integrateForward<decltype(functionImag), Vector>), &functionImag, dim, toArray(xmin), toArray(xmax), 0, prec_rel, prec_abs, ERROR_INDIVIDUAL, &resultImag, &err);

        return resultReal + resultImag * 1i;
      }

      /**
       * Implementation of IntegrationInterface::integrate() function for real valued integrand.
       */
      template <class Function, class Vector, typename std::enable_if_t<!is_complex<typename std::result_of<Function(Vector)>::type>{}>* = nullptr>
      static typename std::result_of<Function(Vector)>::type integrate(Function&& function, Vector xmin, Vector xmax, double prec_rel = 1e-15, double prec_abs = 1e-15, int level = 0) {
      int dim = checkLimits(xmin, xmax);

      if (checkIfIntegrationFirstStepOk(function, xmin, xmax, dim, prec_rel, prec_abs, level) == false) {
        return splitIntegrate(function, xmin, xmax, dim, prec_rel, prec_abs, level+1);
      }

      double result, err;
      pcubature_v(1, reinterpret_cast<cubatureFunctionPointerType>(integrateForward<decltype(function), Vector>), &function, dim, toArray(xmin), toArray(xmax), 0, prec_rel, prec_abs, ERROR_INDIVIDUAL, &result, &err);

      return result;
      }


    private:
      //using cubatureFunctionPointerType = int(*)(unsigned int, long long unsigned int, const double*, void*, unsigned int, double*);
      using cubatureFunctionPointerType = integrand_v;

      // pcubature first evaluates function at xmin, xmax and 0.5(xmin+xmax)
      // if function value at these points is the same it just returns it times volume of (xmin, xmax) area
      // this function checks for such situation
      template <class Function, class Vector>
      static bool checkIfIntegrationFirstStepOk(Function&& function, Vector xmin, Vector xmax, int dim, double prec_rel, double prec_abs, int level) {
        if (level > 5) {
          return true;
        }
        Vector xmid = xmin;
        for (int i = 0; i < dim; ++i) {
          xmid[i] = 0.5 * (xmin[i] + xmax[i]);
        }
        if (std::abs(function(xmin) - function(xmid)) < 10*prec_abs && std::abs(function(xmid) - function(xmax)) < 10*prec_abs) {
          return false;
        } else {
          return true;
        }
      }

      // this function splits integral into 2**dim at point (xmin[i] + (\lambda_i)*(xmax[i]-xmin[i]) )
      // in order to get correct integral using pcubature function (see checkIfIntegrationFirstStepOk())
      // lambda_i are random numbers in [0.4, 0.6] range
      template <class Function, class Vector>
      static auto splitIntegrate(Function&& function, Vector xmin, Vector xmax, int dim, double prec_rel, double prec_abs, int level) {
        std::random_device r;
        std::default_random_engine e1(r());
        std::uniform_real_distribution<double> uniform_dist(0.4, 0.6);

        typename std::result_of<Function(Vector)>::type integral = 0;

        Vector xmid = xmin;
        for (int i = 0; i < dim; ++i) {
          xmid[i] = xmin[i] + uniform_dist(e1) * (xmax[i] - xmin[i]);
        }

        for (int i = 0; i < pow(2, dim); ++i) {
          int ih = i;
          auto xminh = xmin;
          auto xmaxh = xmax;
          for (int j = 0; j < dim; ++j) {
            if (ih % 2  == 0) {
              xmaxh[j] = xmid[j];
            } else {
              xminh[j] = xmid[j];
            }
            ih /= 2;
          }
          integral += integrate(std::ref(function), xminh, xmaxh, prec_rel, prec_abs, level+1);
        }

        return integral;
      }

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
      static int integrateForward(unsigned ndim, long unsigned npts, const double* xArray, void* fdata, unsigned fdim, double* fval);
  };

  template <class Function, class Vector>
  int CubatureIntegration::integrateForward(unsigned ndim, long unsigned npts, const double* xArray, void* fdata, unsigned fdim, double* fval) {
    using Base_Vector = typename std::remove_cv<typename std::remove_reference<Vector>::type>::type;
    using Base_Function = typename std::remove_cv<typename std::remove_reference<Function>::type>::type;

    for (unsigned j = 0; j < npts; ++j) {
      Base_Vector x(ndim);
      for (long unsigned i = 0; i < ndim; ++i) {
        x[i] = xArray[j*ndim+i];
      }
      fval[j] = (*reinterpret_cast<Base_Function*&&>(fdata))(x);
      //fval[j] = (*reinterpret_cast<Function*&&>(fdata))(x);
    }

    return 0;
  }
}

#endif
