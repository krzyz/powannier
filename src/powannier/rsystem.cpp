#include "rsystem.h"
#include "helpers.h"

namespace POWannier {
  RSystem::RSystem(std::shared_ptr<BlochSystem> bs) : 
    N(bs->N),
    dim(bs->dim),
    cutoff(bs->cutoff),
    _bs(std::move(bs)) {
    generateRs();
  }

  void RSystem::generateRs() {
    generateR1();
  }

  void RSystem::generateR1() {
    int rdim = std::pow(N, dim);
    arma::vec r1Eigenvalues(rdim);
    arma::cx_mat r1Eigenvectors(rdim, rdim, arma::fill::zeros);

    loopThroughMN(0,
        [&] (int i, auto&& loopInner) {
          arma::cx_mat r1dPart(N, N);
          for (size_t j = 0; j < N; ++j) {
            if (N%2 == 0) {
              r1dPart(j, j) = 0.5/std::pow(N, dim);
            } else {
              r1dPart(j, j) = 0;
            }
          }

          loopInner([&] (NPoint m, NPoint mp, arma::cx_double el) {
              r1dPart(m(0), mp(0)) = el;
          });

          arma::vec eigval;
          arma::cx_mat eigvec;

          arma::eig_sym(eigval, eigvec, r1dPart);

          r1Eigenvalues.subvec(i*N, (i+1)*N - 1) = std::pow(N, dim) * std::move(eigval);
          r1Eigenvectors.submat(i*N, i*N, (i+1)*N - 1, (i+1)*N - 1) = std::move(eigvec);
        });

    _riEigenvalues.push_back(r1Eigenvalues);
    _riEigenvectors.push_back(r1Eigenvectors);
  }

  void RSystem::generateSubR(arma::cx_mat eigenvalues, int innerDim) {

  }


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

      auto&& loopInner = [&] (auto&& afterNs) {
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

                  double fexp = (ns[nj](0) - ns[ni](0)) * N + (mj - mi);
                  if (fexp != 0) {
                    arma::cx_double phase = std::pow(-1, 1 + fexp);
                    if (N%2 == 0) {
                      phase *= std::exp(std::complex<double>(0, fexp/N*pi));
                    }
                    el += std::conj(_bs->blochC(m, n)) *
                          _bs->blochC(mp, np) * phase *
                          std::complex<double>(0, 1.0/(2*pi*fexp));
                  }
                }
              }
            }

            afterNs(m, mp, el);
          }
        }
      };

      func(i, loopInner);
    }
  }

  Wannier RSystem::getWannier(NPoint n) {
    const arma::vec& r1Eigenvalues = _riEigenvalues[0];
    std::cout << r1Eigenvalues << std::endl;
    int size = r1Eigenvalues.size();
    arma::uvec inds = arma::linspace<arma::uvec>(0, size-1, size);

    NPoint pos = n;
    if (N%2 == 0) {
      pos.for_each([&] (auto& x) {x += N/2 - 1;});
    } else {
      pos.for_each([&] (auto& x) {x += (N-1) / 2;});
    }

    std::cout << pos << std::endl;
    return Wannier(_bs, _riEigenvectors[0].col(pos(0)));
  }
}
