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
  class PotentialBuilder;

  /**
   * @brief Representation of a potential function containing lattice (@f$
   * \vec{a_i} @f$) and reciprocal lattice (@f$ \vec{b_i} @f$) bases vectors.
   * 
   * The periodic potential can be expressed in a so-called plane wave
   * expansion as: @f$ V(\vec{r}) = \sum_n V_n e^{i (\sum_j n_j \vec{b_j})
   * \vec{r}} @f$. @p n is a set of indices @f$ n_i @f$. The $V_n$
   * coefficients are called Fourier coefficient throughout the documentation.
   * 
   */
  class Potential {
    friend class PotentialBuilder;
    friend class PotentialFourierBuilderAddition;
    friend class PotentialFunctionBuilder;

    public:
      /// Default copy constructor.
      Potential(const Potential&) = default;
      /// Default move constructor.
      Potential(Potential&&) = default;

      /**
       * @brief Initialization of the Potential building proces
       * @returns
       *   Builder for a Potential class.
       */
      static PotentialBuilder create();

      /**
       * @brief Value of the potential at a given position.
       */
      double operator()(Position r) const;

      /**
       * @brief Simple getter.
       * @returns
       *   The dimension of the given potential.
       */
      const int dim() const;

      /**
       * @brief Simple getter.
       * @returns
       *   Lattice basis vectors @f$ \vec{a_i} @f$.
       */
      const LatticeBasis& latticeBasis() const;

      /**
       * @brief Simple getter.
       * @returns
       *   Reciprocal lattice basis vectors @f$ \vec{b_i} @f$.
       */
      const ReciprocalBasis& reciprocalBasis() const;

      /**
       * @brief Simple getter.
       * @returns
       *   Fourier coefficients as pairs of @f$ n @f$ and @f$ V_n @f$.
       */
      const FourierCoefficients& fourierCoefficients() const;
      /**
       * @brief Simple getter.
       * @returns
       *  Volume of an elementary cell.
       */
      double elementaryCellVolume() const;

    private:
      int _dim;
      LatticeBasis _latticeBasis;
      ReciprocalBasis _reciprocalBasis;
      FourierCoefficients _coefficients;

      Potential() {}
  };
}

#endif
