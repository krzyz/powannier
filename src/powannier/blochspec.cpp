#include "blochspec.h"

namespace POWannier {
  BlochSpec::BlochSpec(const Potential& V, double kl, int cutoff, double s) :
    kl(kl), 
    cutoff(cutoff),
    dim(V.dim()),
    s(s),
    ns(nspace(cutoff, dim)),
    reciprocalBasis(V.reciprocalBasis()) {}

  std::tuple<arma::vec, arma::cx_mat>
      BlochSpec::generate(const Potential& V, ReciprocalVector k) const {
    arma::cx_mat hamiltonian(ns.size(), ns.size(), arma::fill::zeros);

    for (size_t i = 0; i < ns.size(); ++i) {
      const NPoint& n = ns[i];

      ReciprocalVector nab = ( k + reciprocalBasis * n) * (1.0/kl);

      hamiltonian(i, i) += arma::dot(nab, nab);

      for (const auto& coefficient : V.fourierCoefficients()) {
        NPoint nV = coefficient.first;
        arma::cx_double valueV = coefficient.second;

        NPoint np = n - nV;
        int j = nIndex(np, cutoff);

        if (j >= 0 && j < ns.size()) {
          hamiltonian(i, j) += s * valueV;
        }
      }
    }

    arma::vec energies;
    arma::cx_mat eigenvectors;
    arma::eig_sym(energies, eigenvectors, hamiltonian);

    return std::make_tuple(std::move(energies), std::move(eigenvectors));
  }

  std::complex<double> BlochSpec::blochC(const arma::cx_mat& eigenvectors,
      NPoint n, int band) const {
    return eigenvectors(nIndex(n, cutoff), band);
  }

  std::complex<double> BlochSpec::bloch(const arma::cx_mat& eigenvectors,
      ReciprocalVector k, Position r, int band) const {
    std::complex<double> value = 0;

    auto& coefficients = eigenvectors.col(band);

    for (size_t i = 0; i < ns.size(); ++i) {
      ReciprocalVector g = reciprocalBasis * ns[i];
      value += blochC(eigenvectors, ns[i], band) *
          std::exp(std::complex<double>(0, arma::dot(g, r)));
    }

    value *= std::exp(std::complex<double>(0, arma::dot(k, r)));

    return value;
  }
}