#ifndef POWANNIER_POTENTIALBUILDER_H
#define POWANNIER_POTENTIALBUILDER_H

#include <vector>
#include <queue>
#include <stdexcept>

#include "../potential.h"

namespace POWannier {
  class Potential;
  class PotentialFourierBuilder;
  class PotentialFunctionBuilder;

  /**
   * @brief Main builder for the Potential class.
   * 
   * Basic usage code sample can be found in Potential class description.
   * After providing the basis vectors using setBasis() member function,
   * there are currently two methods of defining the potential:
   * - First of these is simply calling evaluateFromFunction() with 
   *   some callable object taking real-valued <tt>n</tt>-dimensional
   *   vectors and returning real values. In that case, as the end of the
   *   building process, the fourier coefficients will be calculated using
   *   a given IntegrationProvider (which by default is CubatureIntegration).
   *   Several options may be set before the final call using methods of
   *   PotentialFunctionBuilder class.
   * - Second one, enabled by calling addFourierCoefficients, consists of
   *   supplying the coefficients calculated by other means (using the
   *   methods of PotentialFourierBuilder class).
   * 
   * In both of these cases finalization of the building procces is made
   * by a call to complete() function, which returns final Potential object.
   * 
   * @note
   *   The call to setBasis() is mandatory. \n
   *   The only object containing the intermediate result of the build process
   *   (and not only the references to it) is the PotentialBuilder instance 
   *   created by the call to Potential::create(). Because of that, the whole
   *   building proces either has to be finished by a single chain of method 
   *   calls, or the result of the call to Potential::create() has to be stored
   *   in a named variable, e.g.:
   *   @code
   *   // definition of basis, potential funtion, etc.
   * 
   *   // Option 1
   *   auto V = Potential.create()
   *     .setBasis(basis)
   *     .evaluateFromFunction(func)
   *     .complete();
   *   
   *   // Option 2
   *   PotentialBuilder builder = Potential::create();
   *   builder.setBasis(basis);
   *   auto fourierBuilder = builder.addFourierCoefficients();
   *   // calls to fourier Builder
   *   auto V = fourierBuilder.complete();
   *   @endcode
   */
  class PotentialBuilder {
    public:

      /**
       * @brief Set the basis vectors of the potential.
       * @param basis
       *   Basis vectors, stored in std::vector container.
       * @returns
       *   Reference to the current builder.
       */
      PotentialBuilder& setBasis(std::vector<POWannier::Vector> basis);

      /**
       * @brief Start directly adding values of Fourier coefficients of the
       *   plane wave expansion
       * @returns
       *   Builder of type PotentialFourierBuilder.
       */
      PotentialFourierBuilder addFourierCoefficients();

      /**
       * @brief Start the process of calculating the Fourier coefficients
       *   by the means of direct integration of the potential function.
       * @param func
       *   The function of the potential.
       * @returns
       *   Builder of type PotentialFunctionBuilder.
       * @pre
       *   @p func has to take Position as an argument and return real values.
       */
      // the implementation of this template is in potentialfunctionbuilder.h
      template <class Function>
      PotentialFunctionBuilder evaluateFromFunction(Function&& func);
      
      /**
       * @brief Finish building the Potential and return the result.
       * @returns
       *   Finished object of Potential class.
       */
      Potential&& complete();


    protected:
      /// A reference to the Potential object that is currently being built.
      Potential& _potential;

      /**
       * @brief Constructor taking a reference to Potential.
       * 
       * Shallow copy of the given Potential - only the reference is stored.
       * Not to be used to create a new builder.
       * @param potential
       *   Potential to use in initializing the builder.
       */
      explicit PotentialBuilder(Potential& potential) :
        _potential(potential) {}

    private:
      friend class Potential;

      enum class State {
        BasisNotDeclared,
        BasisDeclared,
      };

      PotentialBuilder();

      void ensureBasisDefined();

      Potential _p;
      State _state;
  };

 }

#endif
