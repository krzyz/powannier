#ifndef _POWANNIER_POTENTIALFOURIERBUILDERADDITION_H
#define _POWANNIER_POTENTIALFOURIERBUILDERADDITION_H

#include "potentialbuilder.h"

namespace POWannier {
  /**
   * @brief Builder of the Potential class directly from Fourier coefficients - 
   *   waiting for values state.
   */
  class PotentialFourierBuilderAddition : private PotentialBuilder {
    public:
      /**
       * @brief Provide batch of Fourier coefficients values for indices
       * defined earlier.
       * @param values
       *   List of values corresponding to indices supplied earlier (in order).
       * @returns
       *   Builder of type PotentialFourierBuilder.
       */
      PotentialFourierBuilder values(std::vector<Complex> values);

      /**
       * @brief Provide a single Fourier coefficient value for the index
       * defined earlier.
       * @param value
       *  Value corresponding to the index supplied earlier (next in order).
       * @returns
       *   Builder of type PotentialFourierBuilder.
       */
      PotentialFourierBuilder value(Complex value);

    private:
      friend class PotentialFourierBuilder;
      explicit PotentialFourierBuilderAddition(
        Potential& potential, std::queue<NPoint> indices) :
        PotentialBuilder(potential), _indices(indices) {}
      
      void addToPotentialCoefficients(Complex value);
      std::queue<NPoint> _indices;
  };
}

#endif