#include "blochsystem.h"
#include "helpers.h"


namespace POWannier {
  BlochSystem::BlochSystem(std::shared_ptr<const Potential> V, double kl, int cutoff, int N, double s) :
    kl(kl), 
    cutoff(cutoff),
    dim(V->dim),
    N(N),
    s(s),
    ms(mspace(N, dim)),
    ns(nspace(cutoff, dim)),
    V(std::move(V))
     {
      _reciprocalBasis = this->V->reciprocalBasis();
      generateAll();
    }

  void BlochSystem::generateAll() {
    _energies = std::vector<arma::vec>(ms.size());
    _eigenvectors = std::vector<arma::cx_mat>(ms.size());

    for (size_t mi = 0; mi < ms.size(); ++mi) {
      const auto& m = ms[mi];
      ReciprocalVector k = kFromM(m);

      arma::cx_mat hamiltonian(ns.size(), ns.size(), arma::fill::zeros);

      for (size_t i = 0; i < ns.size(); ++i) {
        NPoint& n = ns[i];

        ReciprocalVector nab = ( k + _reciprocalBasis * n) * (1.0/kl);

        hamiltonian(i, i) += arma::dot(nab, nab);

        for (const auto& coefficient : V->fourierCoefficients()) {
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
      /*
      for (auto mx : m) {
        std::cout << mx << " ";
      }

      for (int i = 0; i < 16; ++i) {
        std::cout << energies[i] << " ";
      }
      std::cout << std::endl;
      */

      _energies[mi] = energies;
      _eigenvectors[mi] = eigenvectors;
    }
  }

  const arma::vec& BlochSystem::energies(NPoint m) const {
    return _energies[mIndex(m, N)];
  }
  
  const arma::cx_mat& BlochSystem::eigenvectors(NPoint m) const {
    return _eigenvectors[mIndex(m, N)];
  }

  std::complex<double> BlochSystem::blochC(NPoint m, NPoint n, int band) const {
    return eigenvectors(m)(nIndex(n, cutoff), band);
  }

  std::complex<double> BlochSystem::bloch(NPoint m, Position r, int band) const {
    ReciprocalVector k = kFromM(m);
    std::complex<double> value = 0;

    auto& coefficients = eigenvectors(m).col(band);

    for (size_t i = 0; i < ns.size(); ++i) {
      ReciprocalVector g = _reciprocalBasis * ns[i];
      value += blochC(m, ns[i], band) * std::exp(std::complex<double>(0, arma::dot(g, r)));
    }

    value *= std::exp(std::complex<double>(0, arma::dot(k, r)));

    return value;
  }


  ReciprocalVector BlochSystem::kFromM(NPoint m) const {
    return _reciprocalBasis * m / N;
  }
}
