
namespace POWannier {
  const double pi = arma::datum::pi;
  using Position = arma::rowvec;
  using Vector = arma::rowvec;
  using ReciprocalVector = arma::vec;
  using ReciprocalBasis = arma::mat;
  using LatticeBasis = arma::mat;
  using NPoint = arma::ivec;
  using FourierCoefficients = std::vector<std::pair<NPoint, arma::cx_double>>;
}
