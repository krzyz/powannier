#ifndef _POWANNIER_POTENTIALFOURIERBUILDER_H
#define _POWANNIER_POTENTIALFOURIERBUILDER_H

#include "potentialbuilder.h"

namespace POWannier {
  class PotentialFourierBuilderAddition;
  /**
   * @brief Builder of the Potential class directly from Fourier coefficients - base state.
   */
  class PotentialFourierBuilder : private PotentialBuilder {
    public:
      /**
       * @brief Begin adding Fourier coefficients by specyfing
       *   which indices will be provided next.
       * @param indices
       *   List of <tt>n</tt>-dimensional indices.
       * @returns
       *   Builder of type PotentialFourierBuilderAddition.
       */
      PotentialFourierBuilderAddition provideMultiple(std::vector<NPoint> indices);

      /**
       * @brief Begin adding Fourier coefficients by specyfing
       *   which index will be provided next.
       * @param index
       *   <tt>n</tt>-dimensional index.
       * @returns
       *   Builder of type PotentialFourierBuilderAddition.
       */
      PotentialFourierBuilderAddition provideOne(NPoint index);

      /**
       * @brief Finish building the Potential and return the result.
       * @returns
       *   Finished object of Potential class.
       */
      using PotentialBuilder::complete;
    private:
      friend class PotentialBuilder;
      friend class PotentialFourierBuilderAddition;
      explicit PotentialFourierBuilder(Potential& potential) :
        PotentialBuilder(potential) {
      }
  };
}

#endif