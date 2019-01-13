#ifndef POWANNIER_ALIASES_H
#define POWANNIER_ALIASES_H

#include <armadillo>

namespace POWannier {
  const double pi = arma::datum::pi;
  using Complex = arma::cx_double;
  using Position = arma::rowvec;
  using Vector = arma::rowvec;
  using ReciprocalVector = arma::vec;
  using ReciprocalBasis = arma::mat;
  using LatticeBasis = arma::mat;
  using NPoint = arma::ivec;
  using FourierCoefficients = std::vector<std::pair<NPoint, Complex>>;
  using Eigensystem = std::tuple<arma::vec, arma::cx_mat>;
}

#endif