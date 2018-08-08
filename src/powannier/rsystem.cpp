#include "rsystem.h"
#include "helpers.h"

namespace POWannier {
  RSystem::RSystem(std::shared_ptr<BlochSystem> bs) : 
    N(bs->N),
    dim(bs->dim),
    cutoff(bs->cutoff),
    _bs(std::move(bs)) {
    _r1EigenSystem = generateR1();
  }

/*
  void RSystem::test() {
    int rdim = std::pow(N, dim);
    arma::cx_mat r(rdim, rdim, arma::fill::zeros);
    loopThroughMN(0, auto&& loopInner[&] (NPoint m, NPoint mp, arma::cx_double el) {
      int mi = mIndex(m, N);
      int mj = mIndex(mp, N);
      r(mi, mj) = arma::cx_double(mi, mj);
    });
    r.print();

    loopThroughMN(1, [&] (NPoint m, NPoint mp, arma::cx_double el) {
      int mi = mIndex(m, N);
      int mj = mIndex(mp, N);
      r(mi, mj) = arma::cx_double(mi, mj);
    });
    r.print();

  }
  */

  Wannier RSystem::getWannier(NPoint n) {
    //test();
    NPoint pos = n;
    if (N%2 == 0) {
      pos.for_each([&] (auto& x) {x += N/2 - 1;});
    } else {
      pos.for_each([&] (auto& x) {x += (N-1) / 2;});
    }

    arma::cx_mat eigenvectors = r1Eigenvectors();
    arma::vec eigenvalues = r1Eigenvalues();

    int idim = std::pow(N, dim);

    arma::cx_mat subspace(idim, idim, arma::fill::eye);

    for (int i = 0; i < dim-1; ++i) {
      idim = std::pow(N, dim-1-i);
      std::cout << "idim: " << idim << std::endl;

      subspace = subspace * eigenvectors.cols(idim*pos(i), idim * (pos(i)+1) - 1);
      std::cout << "subspace eigenvalues: " << std::endl << eigenvalues.subvec(idim*pos(i), idim * (pos(i)+1) - 1);

      EigenSystem subsystem = generateSubR(subspace, i+1);

      eigenvalues = std::move(std::get<0>(subsystem));

      std::cout << "r2 eigenvalues: " << std::endl << eigenvalues;
      eigenvectors = std::move(std::get<1>(subsystem));
    }

    arma::cx_mat coefficients = subspace * eigenvectors.col(pos(dim-1));

    return Wannier(_bs, coefficients);
  }

  const arma::cx_mat& RSystem::r1Eigenvectors() {
    return std::get<1>(_r1EigenSystem);
  }

  const arma::vec& RSystem::r1Eigenvalues() {
    return std::get<0>(_r1EigenSystem);
  }

  void RSystem::sortEigensystem(arma::vec& eigenvalues, arma::cx_mat& eigenvectors) {
    arma::uvec sortedIndices = arma::sort_index(eigenvalues);
    eigenvalues = eigenvalues(sortedIndices);
    eigenvectors = eigenvectors.cols(sortedIndices);
  }

  arma::uvec RSystem::rearrangedIndicesR1() {
    int lddim = std::pow(N, dim-1);
    arma::uvec indices(std::pow(N, dim));
    arma::uvec oneMIndices = arma::linspace<arma::uvec>(0, std::pow(N, dim)-N, lddim);
    for (int i = 0; i < N; ++i) {
      indices.subvec(i*lddim, (i+1)*lddim-1) = i+oneMIndices;
    }

    std::cout << oneMIndices << std::endl;
    std::cout << indices << std::endl;
    return indices;
  }

  EigenSystem RSystem::generateR1() {
    int rdim = std::pow(N, dim);
    arma::vec r1Eigenvalues(rdim);
    arma::cx_mat r1Eigenvectors(rdim, rdim, arma::fill::zeros);

    loopThroughMN(0,
      [&] (int i, auto&& loopInner) {
        arma::cx_mat r1dPart(N, N);
        /*
        for (int j = 0; j < N; ++j) {
          if (N%2 == 0) {
            r1dPart(j, j) = 0.5;
          } else {
            r1dPart(j, j) = 0;
          }
        }
        */

        loopInner([&] (NPoint m, NPoint mp, arma::cx_double el) {
            r1dPart(m(0), mp(0)) = el;
        });

        arma::vec eigval;
        arma::cx_mat eigvec;

        arma::eig_sym(eigval, eigvec, r1dPart);

        r1Eigenvalues.subvec(i*N, (i+1)*N - 1) = std::move(eigval);
        r1Eigenvectors.submat(i*N, i*N, (i+1)*N - 1, (i+1)*N - 1) = std::move(eigvec);
      });

    sortEigensystem(r1Eigenvalues, r1Eigenvectors);
    auto rearrangedIndices = rearrangedIndicesR1();
    r1Eigenvectors = r1Eigenvectors.rows(rearrangedIndices);

    std::cout << "r1 eigenvalues:" << std::endl << r1Eigenvalues << std::endl;

    return std::make_tuple(std::move(r1Eigenvalues), std::move(r1Eigenvectors));
  }

  EigenSystem RSystem::generateSubR(arma::cx_mat eigenvectors, int innerDim) {
    int rdim = eigenvectors.n_cols;
    arma::cx_mat subr(rdim, rdim, arma::fill::zeros);

/*
    for (size_t mi = 0; mi < _bs->ms.size(); ++mi) {
      for (int ei = 0; ei < rdim; ++ei) {
        for (int ej = ei; ej < rdim; ++ej) {
          subr(ei, ej) += 0.5 * std::conj(eigenvectors(mi, ei)) *
                          eigenvectors(mi, ej);
        }
      }
    }
    */

    loopThroughMN(innerDim, 
      [&] (int i, auto&& loopInner) {
        loopInner([&] (NPoint m, NPoint mp, arma::cx_double el) {
          int mi = mIndex(m, N);
          int mj = mIndex(mp, N);
          std::cout << mi << " " << mj << std::endl;
          for (int ei = 0; ei < rdim; ++ei) {
            for (int ej = ei; ej < rdim; ++ej) {
              subr(ei, ej) += el * std::conj(eigenvectors(mi, ei)) *
                            eigenvectors(mj, ej);
              if (mi != mj) {
                subr(ei, ej) += std::conj(el * eigenvectors(mj, ei)) *
                              eigenvectors(mi, ej);
              }
            }
          }
        });
      });

    arma::vec eigval;
    arma::cx_mat eigvec;

    arma::eig_sym(eigval, eigvec, subr);

    return make_tuple(std::move(eigval), std::move(eigvec));
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
          m(innerDim) = mi;

          if (N%2 == 0) {
            afterNs(m, m, 0.5);
          } else {
            afterNs(m, m, 0.0);
          }

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

            afterNs(m, mp, static_cast<double>(N) * el);
          }
        }
      };

      func(i, loopInner);
    }
  }

}
