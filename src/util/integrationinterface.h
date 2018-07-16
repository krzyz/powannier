#ifndef INTEGRATIONINTERFACE
#define INTEGRATIONINTERFACE


namespace POWannier {

  /**
   * \brief An interface providing a static function to calculate integrals.
   * Classes implementing this interface have a static function integrate which can be used to calculate one or multidimensional definite integrals.
   */

  class IntegrationInterface {
    public:
      virtual ~IntegrationInterface() {}

      /**
       * Calculate a definite integral of a provided function \f$ f \f$ over a given period/area.
       *
       * @param function A function to integrate. Can be a lambda expression, functor or a pointer to a function.
       * @param xmin A vector containing lower limit of the integral.
       * @param xmax A vector containing upper limit of the integral.
       *
       * @pre @p xmin and @p xmax must have the same size \f$ n \f$.
       * @pre @p function must return real or complex number for a given \f$ n \f$ dimensional vector.
       *
       * @returns The result of the integration: \f$ \int dx^n f(x_1, \ldots, x_n) \f$.
       */
      template <class Function, class Vector>
      static auto integrate(
          Function&& function, Vector&& xmin, Vector&& xmax 
        ) -> decltype(function(xmin));
  };
}

#endif
