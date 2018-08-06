#include "rsystem.h"
#include "helpers.h"

namespace POWannier {
  RSystem::RSystem(std::shared_ptr<BlochSystem> bs) : 
    N(bs->N),
    dim(bs->dim),
    cutoff(bs->cutoff),
    _bs(std::move(bs)) {
    generateR();
  }

  void RSystem::generateR() {
    /*
    auto ns = nspace(cutoff, 1);
    int bdim = ns.size();
    std::vector<int> bands = {0, 1};
    //int dim = bands.size() * std::pow(bdim, dim) ;
    int rdim = std::pow(N, dim);

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
              n(0) = ns[ni](0);
              for (int nj = 0; nj < bdim; ++nj) {
                np(0) = ns[nj](0);

                double fexp = ns[nj](0) - ns[ni](0) + (mj - mi) * N;
                if (fexp != 0) {
                  el += std::conj(_bs->blochC(m, n)) * _bs->blochC(m, np) *
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
      r1Eigenvectors.submat( i*N, i*N, (i+1) * N - 1, (i+1) * N - 1) = std::move(eigvec);
    }
    */

    int rdim = std::pow(N, dim);
    arma::vec r1Eigenvalues(rdim);
    arma::cx_mat r1Eigenvectors(rdim, rdim, arma::fill::zeros);

    loopThroughMN(0,
        [&] (int i, auto&& loopInner) {
          arma::cx_mat r1dPart(N, N);
          for (size_t j = 0; j < N; ++j) {
          r1dPart(j, j) = 0;
          }

          loopInner([&] (NPoint m, NPoint mp, arma::cx_double el) {
              r1dPart(m(0), mp(0)) = el;
          });

          arma::vec eigval;
          arma::cx_mat eigvec;
          arma::eig_sym(eigval, eigvec, r1dPart);

          r1Eigenvalues.subvec(i * N, (i+1) * N -1) = std::pow(N, dim) * std::move(eigval);
          r1Eigenvectors.submat( i*N, i*N, (i+1) * N - 1, (i+1) * N - 1) = std::move(eigvec);
        });

    std::cout << r1Eigenvalues << std::endl;
    std::cout << r1Eigenvectors << std::endl;
  }

  /*
  void RSystem::generateSubR(arma::cx_mat eigenvalues, int innerDim) {

  }
  */


  void RSystem::loopThroughMN(int innerDim, auto&& func) {

    auto ns = nspace(cutoff, 1);
    int bdim = ns.size();
    std::vector<int> bands = {0, 1};

    auto outerMSpace = mspace(N, dim-1);
    auto outerNSpace = nspace(cutoff, dim-1);

    arma::uvec outerDims = arma::linspace<arma::uvec>(0, dim-1, dim);
    outerDims.shed_row(innerDim);

    for (int i = 0; i < outerMSpace.size(); ++i) {

      NPoint m(dim);
      m(outerDims) = outerMSpace[i];
      NPoint mp = m;

      func(i, [&] (auto&& afterNs) {
          for (int mi = 0; mi < N; ++mi) {
            m(innerDim) = mi ;
            for (int mj = mi+1; mj < N; ++mj) {
              mp(innerDim) = mj;
              arma::cx_double el;

              for (int j = 0; j < outerNSpace.size(); ++j) {
                NPoint n(dim);
                n(outerDims) = outerNSpace[j];
                NPoint np = n;

                for (int ni = 0; ni < bdim; ++ni) {
                  n(innerDim) = ns[ni](0);
                  for (int nj = 0; nj < bdim; ++nj) {
                    np(innerDim) = ns[nj](0);

                    double fexp = ns[nj](0) - ns[ni](0) + (mj - mi) * N;
                    if (fexp != 0) {
                      el += std::conj(_bs->blochC(m, n)) *
                            _bs->blochC(m, np) *
                            std::complex<double>(0, -1.0/fexp);
                    }
                  }
                }
              }

              afterNs(m, mp, el);
            }
          }
      });
    }
  }

}
