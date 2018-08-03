#include "rsystem.h"
#include "helpers.h"

namespace POWannier {
  void RSystem::generateRHamiltonians(BlochSystem bs) {
    int N = 3;
    int dim = 2;
    int cutoff = 1;
    auto ns= nspace(cutoff, 1);
    int bdim = ns.size();
    std::vector<int> bands = {0, 1};
    //int dim = bands.size() * std::pow(bdim, dim) ;
    int rdim = std::pow(bdim, dim);

    arma::vec r1Eigenvalues(rdim);
    arma::cx_mat r1Eigenvectors(rdim, rdim, arma::fill::zeros);

    auto outerMSpace = mspace(N, dim-1);
    auto outerNSpace = nspace(cutoff, dim-1);
    for (int i = 0; i < outerMSpace.size(); ++i) {
      arma::cx_mat r1dPart(N, N);
      NPoint m(dim);
      m.tail(dim-1) = outerMSpace[i];
      NPoint mp = m;

      for (int mi = 0; mi < N; ++mi) {
        m(0) = mi ;
        for (int mj = mi+1; mj < N; ++mj) {
          mp(0) = mj;
          arma::cx_double el;

          for (int j = 0; j < outerNSpace.size(); ++j) {
            NPoint n(dim);
            n.tail(dim-1) = outerNSpace[j];
            NPoint np = n;

            for (int ni = 0; ni < bdim; ++ni) {
              n(0) = ni;
              for (int nj = 0; nj < bdim; ++nj) {
                np(0) = nj;

                double fexp = ns[nj](0) - ns[ni](0);
                if (fexp != 0) {
                  el += std::conj(bs.blochC(m, n)) * bs.blochC(m, np) *
                        std::complex<double>(0, -1.0/fexp);
                }
              }
            }
          }

          r1dPart(mi, mj) = el;
        }
      }


      arma::vec eigval;
      arma::cx_mat eigvec;
      arma::eig_sym(eigval, eigvec, r1dPart);

      r1Eigenvalues.subvec(i * N, (i+1) * N -1) = std::move(eigval);
      r1Eigenvectors.submat( i*N, (i+1) * N - 1, i*N, (i+1) * N - 1) = std::move(eigvec);
    }

    std::cout << r1Eigenvalues << std::endl;
    std::cout << r1Eigenvectors << std::endl;
  }
}