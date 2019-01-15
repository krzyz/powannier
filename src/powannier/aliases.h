#ifndef POWANNIER_ALIASES_H
#define POWANNIER_ALIASES_H

#include <armadillo>

/**
 * @brief
 *   The default namespace of the project, which in addition to all the classes
 *   contains typedefs, enums and helper functions.
 */
namespace POWannier {
  /// Pi constant (approx. 3.14).
  const double pi = arma::datum::pi;
  /// Complex number type
  using Complex = arma::cx_double;
  /**
   * @brief position type (both in real space and in lattice basis 
   * vectors space), row vector.
   */
  using Position = arma::rowvec;

  /**
   * @brief
   *   Similar to Position (however here it is used for values 
       not related to system's set of coordinates), row vector.
   */
  using Vector = arma::rowvec;
  /// vector or position in reciprocal space, column vector.
  using ReciprocalVector = arma::vec;
  /// Set of reciprocal basis vectors (the vectors are columns here).
  using ReciprocalBasis = arma::mat;
  /// Set of lattice basis vectors (the vectors are rows here).
  using LatticeBasis = arma::mat;
  /// Column vector containing a set of integer indices
  using NPoint = arma::ivec;

  /**
   * @brief
   *   Type for storing Fourier series coefficients (used in plane wave
   *   approximations of periodic functions, such as potential function \f$V\f$).
   */
  using FourierCoefficients = std::vector<std::pair<NPoint, Complex>>;

  /**
   * @brief
   *   Tuple for storing both the eigenvalues and the eigenvectors of
   *   Hermitian eigenvalue problems.
   */
  using Eigensystem = std::tuple<arma::vec, arma::cx_mat>;
}

#endif