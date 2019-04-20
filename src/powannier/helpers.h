#ifndef POWANNIER_HELPERS_H
#define POWANNIER_HELPERS_H

#include <complex>
#include <vector>

#include "aliases.h"

namespace POWannier {
  #pragma omp declare reduction \
    (compadd:arma::cx_double:omp_out=omp_out+omp_in) \
    initializer(omp_priv=0)

  /**
   * @brief 
   *   Convenience function returning 0 for real values
   *   smaller than given precision.
   * 
   * @param number
   *   Number to approximate.
   * @param prec
   *   Precision to be used.
   * 
   * @returns
   *   0 if @p |number| < @p prec, \n
   *   @p number otherwise.
   */
  double chop(double number, double prec = 1e-15);

  /**
   * @brief 
   *   Convenience function returning 0 for complex 
   *   values smaller than given precision.
   * 
   * @param number
   *   Number to approximate.
   * @param prec
   *   Precision to be used.
   * 
   * @returns
   *   @p chop(re(number)) + \f$i\f$ @p chop(im(number)).
   */

  std::complex<double> chop(std::complex<double> number, double prec = 1e-15);

  /**
   * @brief
   *   The space of sets of @p dim indices, with each index ranging
   *   from @p -cutoff to @p cutoff.
   * 
   * @param cutoff
   *   The cutoff parameter.
   * @param dim
   *   The dimension of the index set.
   * 
   * @returns 
   *   Vector containing all indices in the range.
   */
  std::vector<NPoint> nspace(int cutoff, int dim);

  /**
   * @brief
   *   The space of sets of @p dim indices, with each index ranging
   *   from @p 0 to @p N-1
   * 
   * @param N 
   *   The range of indices (number of single index values).
   * @param dim
   *   The dimension of the index set.
   * 
   * @returns 
   *   Vector containing all indices in the range.
   */

  std::vector<NPoint> mspace(int N, int dim);

  /**
   * The position of a given index in the @ref nspace of indices.
   * 
   * @param n
   *   Given index.
   * @param cutoff
   *   Cutoff used in @ref nspace.
   * @returns
   *  The position of index @p n in @ref nspace.
   */
  int nIndex(const NPoint& n, int cutoff);

  /**
   * The position of a given index in the @ref mspace of indices.
   * @param n
   *   Given index.
   * @param N
   *   N value used in @ref mspace.
   * @returns
   *  The position of index @p n in @ref mspace.
   */
  int mIndex(const NPoint& n, int N);

  /**
   * @brief
   *   The position of a given index in the basis of indices used
   *   for storing Wannier function coefficients
   * @param m
   *   Given index (of reciprocal space vectors, \f$ \vec{k} = \sum_i (2\pi 
   *   m_i / N) \vec{i} \f$, where \f$ \vec{i} \f$ are the Carthesian coordinate
   *   system basis vectors.
   * @param N
   *   System size (number of cells in each direction).
   * @param band
   *   The band index.
   * @returns
   *   The position of index @p m in a basis used in Wannier
   *   function storage.
   */
  int bmIndex(const NPoint& m, int N, int band);

  std::vector<int> commonBands(Wannier& wan1, Wannier& wan2);
}

#endif