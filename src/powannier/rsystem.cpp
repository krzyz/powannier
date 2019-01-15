#include <stdexcept>
#include "rsystem.h"
#include "helpers.h"

using namespace std::placeholders;

namespace POWannier {
  RSystem::RSystem(std::shared_ptr<BlochSystem> bs, std::vector<int> bands) : 
    N(bs->N),
    dim(bs->dim()),
    mdim(std::pow(N, dim)),
    bmdim(std::pow(N, dim) * bands.size()),
    bands(bands),
    cutoff(bs->cutoff()),
    _bs(std::move(bs)),
    _eigensystemCache(std::bind(&RSystem::calculateEigensystem, this, _1)),
    _subEigensystemCache(std::bind(static_cast<Eigensystem(RSystem::*)(Ipp)>(&RSystem::calculateSubEigensystem), this, _1)),
    _subspaceCache(std::bind(&RSystem::calculateSubspace, this, _1)) {
    std::vector<int> positions(bands.size(), 1);
    _wannierPositions = WannierPositions(positions);
  }

  RSystem::RSystem(std::shared_ptr<BlochSystem> bs, std::vector<int> bands, WannierPositions positions) : 
    RSystem(bs, bands) {
    if (_wannierPositions.descendantLeavesNumber() != static_cast<int>(bands.size())) {
      throw std::runtime_error("Number of specified positions of wannier functions must be equal to bands number!");
    }
    _wannierPositions = positions;
  }

  Wannier RSystem::getWannier(NPoint n) {
    NPoint elCellPositions(dim, arma::fill::zeros);
    return getWannier(n, elCellPositions);
  }

  Wannier RSystem::getWannier(NPoint cellLocation, NPoint wanLocation) {
    auto subspace = getSubspace(dim-1, cellLocation, wanLocation);
    auto subEigenvectors = getSubEigenvectors(dim-1, cellLocation, wanLocation);

    int position = wannierPosition(dim-1, cellLocation, wanLocation);

    arma::cx_mat coefficients = subspace * subEigenvectors.col(position);

    return Wannier(_bs, coefficients, bands);
  }

  int RSystem::wannierPosition(int n, NPoint cellLocation, NPoint wanLocation) {
    return std::get<0>(getSubspaceRange(n, cellLocation, wanLocation, true));
  }

  std::tuple<int,int> RSystem::getSubspaceRange(int n, NPoint cellLocation, NPoint wanLocation, bool last) {
    int idim = std::pow(N, dim-n-1);

    auto pos = cellLocation;

    if (N%2 == 0) {
      pos.for_each([&] (auto& x) {x += N/2 - 1;});
    } else {
      pos.for_each([&] (auto& x) {x += (N-1) / 2;});
    }

    auto currentPositions = _wannierPositions;
    for (int i = 0; i < n; ++i) {
      currentPositions = currentPositions
        .getChild(wanLocation(i));
    }

    int mv = bands.size();
    if (last == true) {
      mv = currentPositions.descendantLeavesNumber();
    }

    int subspaceStart = idim * (mv * pos(n)
      + currentPositions.descendantLeavesNumberLeftTo(wanLocation(n)));

    int subspaceEnd = subspaceStart +
        idim * currentPositions
          .getChild(wanLocation(n))
          .descendantLeavesNumber() - 1;

    return std::make_tuple(subspaceStart, subspaceEnd);
  }

  const Eigensystem& RSystem::getEigensystem(int n) {
    return _eigensystemCache.get(n);
  }


  const Eigensystem& RSystem::getSubEigensystem(int n, NPoint cellLocation, NPoint wanLocation) {
    if (n == 0) {
      return getEigensystem(0);
    }

    Ipp ipp = std::make_tuple(
      n, 
      arma::conv_to<std::vector<int>>::from(cellLocation.head(n)), 
      arma::conv_to<std::vector<int>>::from(wanLocation.head(n))
    );

    return _subEigensystemCache.get(ipp);
  }

  const arma::cx_mat& RSystem::getSubspace(int n, NPoint cellLocation, NPoint wanLocation) {

    Ipp ipp = std::make_tuple(
      n, 
      arma::conv_to<std::vector<int>>::from(cellLocation.head(n)), 
      arma::conv_to<std::vector<int>>::from(wanLocation.head(n))
    );

    return _subspaceCache.get(ipp);
  }

  Eigensystem RSystem::calculateSubEigensystem(Ipp ipp) {
    auto n = std::get<0>(ipp);
    auto cellLocation = arma::conv_to<NPoint>::from(std::get<1>(ipp));
    auto wanLocation = arma::conv_to<NPoint>::from(std::get<2>(ipp));


    const auto& subspace = getSubspace(n, cellLocation, wanLocation);
    return calculateSubEigensystem(n, subspace);
  }

  arma::cx_mat RSystem::calculateSubspace(Ipp ipp) {
    auto n = std::get<0>(ipp);
    auto cellLocation = arma::conv_to<NPoint>::from(std::get<1>(ipp));
    auto wanLocation = arma::conv_to<NPoint>::from(std::get<2>(ipp));

    if (n == 0) {
      return arma::cx_mat(bmdim, bmdim, arma::fill::eye);
    }

    const auto& subspace = getSubspace(n-1, cellLocation, wanLocation);
    const auto& eigenvectors = getSubEigenvectors(n-1, cellLocation, wanLocation);

    auto subspaceRange = getSubspaceRange(n-1,  cellLocation, wanLocation);
    auto subspaceStart = std::get<0>(subspaceRange);
    auto subspaceEnd = std::get<1>(subspaceRange);

    return subspace * eigenvectors.cols(subspaceStart, subspaceEnd);
  }

  // used to calculate permutations transforming bases used in calculating
  // position operator of dimension "inner" to the main basis
  arma::uvec RSystem::transformFromInner(int inner) const {
    // in order to get the correct permutation, for each set
    // of indices m (in the order in which they appear in Eigensystem
    // calculations related to the "inner" dimension) to the indices m,
    // which correspond to the same set of indices in the order in which
    // they appear in the main basis
    if (inner < 0 || inner > dim-1) {
      throw std::runtime_error("Inner dimension index no in range of dimension indices!");
    }

    // the order in the last dimension is the order of the main basis
    if (inner == dim-1) {
      return arma::linspace<arma::uvec>(0, mdim-1, mdim);
    }

    arma::uvec indices(mdim);

    for (int i = 0; i < mdim; ++i) {
      auto m = _bs->ms[i];
      auto mRearranged = m;

      // the m indices have the form:
      // (m_0, m_1, ... m_{inner-1}, m_{inner+1}, ... m_{dim-1}, m_{inner})
      // while mRearranged are:
      // (m_0, m_1, ... m_{inner-1}, m_{inner},   ... m_{dim-2}, m_{dim-1})
      mRearranged.subvec(inner+1, dim-1) = m.subvec(inner, dim-2);
      mRearranged(inner) = m(dim-1);

      auto inew = mIndex(mRearranged, N);
      indices(inew) = i;
    }

    return indices;
  }

  arma::uvec RSystem::transformToBm() const {
    int odim = std::pow(N, dim-1);
    arma::uvec indices(bmdim);

    int i = 0;

    for (int oi = 0; oi < odim; ++oi) {
      for (size_t bi = 0; bi < bands.size(); ++bi) {
        for (int ii = 0; ii < N; ++ii) {
          indices(bi * mdim + oi * N + ii) = i;
          ++i;
        }
      }
    }

    return indices;
  }

  // same as above, but including the bands in the aggregate basis index
  arma::uvec RSystem::transformFromInnerToBm(int inner) const {
    arma::uvec indices = transformToBm();
    arma::uvec indicesFromInner = transformFromInner(inner);
    for (std::size_t i = 0; i < bands.size(); ++i) {
      arma::uvec indicesInOneBand = indices.subvec(i*mdim, (i+1) * mdim - 1);
      indices.subvec(i*mdim, (i+1) * mdim - 1) = indicesInOneBand(indicesFromInner);
    }

    return indices;
  }

  void RSystem::sortEigensystem(arma::vec& eigenvalues, arma::cx_mat& eigenvectors) const {
    arma::uvec sortedIndices = arma::sort_index(eigenvalues);
    eigenvalues = eigenvalues(sortedIndices);
    eigenvectors = eigenvectors.cols(sortedIndices);
  }

  Eigensystem RSystem::calculateEigensystem(int inner) const {
    int submdim = N*bands.size();

    auto outerMSpace = mspace(N, dim-1);
    auto omdim = outerMSpace.size();

    arma::uvec outerDims = arma::linspace<arma::uvec>(0, dim-1, dim);
    outerDims.shed_row(inner);

    arma::vec reigval(bmdim);
    arma::cx_mat reigvec(bmdim, bmdim, arma::fill::zeros);

    for (std::size_t i = 0; i < omdim; ++i) {

      NPoint m(dim);
      m(outerDims) = outerMSpace[i];

      arma::cx_mat r1dBands(submdim, submdim);

      for (std::size_t bi = 0; bi < bands.size(); ++bi) {
        int first_row = bi*N;
        for (std::size_t bj = bi; bj < bands.size(); ++bj) {
          int first_col = bj*N;
          r1dBands.submat(first_row, first_col, arma::size(N, N)) = 
            rInnerMatrix(inner, m, bi, bj);
        }
      }

      arma::vec eigval;
      arma::cx_mat eigvec;

      arma::eig_sym(eigval, eigvec, r1dBands);

      reigval.subvec(i*submdim, arma::size(eigval)) = std::move(eigval);
      reigvec.submat(i*submdim, i*submdim, arma::size(eigvec)) = std::move(eigvec);
    }

    auto indicesInStandardBasis = transformFromInnerToBm(inner);

    reigvec = reigvec.rows(indicesInStandardBasis);

    sortEigensystem(reigval, reigvec);

    return make_tuple(std::move(reigval), std::move(reigvec));
  }

  Eigensystem RSystem::calculateSubEigensystem(int inner, const arma::cx_mat& subspace) {
    const auto& rEigval = getEigenvalues(inner);
    const auto& rEigvec = getEigenvectors(inner);
    arma::cx_mat subr = 
      subspace.t() * 
      rEigvec * arma::diagmat(rEigval) * rEigvec.t() *
      subspace;

    arma::vec eigval;
    arma::cx_mat eigvec;
    arma::eig_sym(eigval, eigvec, subr);

    return std::make_tuple(std::move(eigval), std::move(eigvec));
  }

  // calculate part of the block in the position operator of dimension "inner"
  //   block corresponds to m indices (m_{i} for i = inner)
  //   part corresponds to bands bi (row) and bj (column)
  // m are integers that reproduce wavewectors from first brillouin zone:
  //   k = (2 * \pi * (m - f(N)) / N, where f(n) = 
  arma::cx_mat RSystem::rInnerMatrix(int inner, NPoint m, int bi, int bj) const {
    auto mp = m;

    arma::cx_mat matrix(N, N);
    for (int mi = 0; mi < N; ++mi) {
      m(inner) = mi;

      int mjstart = 0;
      if (bi == bj) {
        matrix(mi, mi) = (N%2 == 0 ? 0.5 : 0);
        mjstart = mi + 1;
      }
      for (int mj = mjstart; mj < N; ++mj) {
        mp(inner) = mj;
        arma::cx_double el;

        matrix(mi, mj) = static_cast<double>(N) * sumOverNs(inner, m, mp, bi, bj);
      }

    }

    return matrix;
  }

  // summation over all fourier coefficents for indices m (row) and mp (column),
  // and bands bands[bi] (row) and bands[bj] (column)
  arma::cx_double RSystem::sumOverNs(int inner, NPoint m, NPoint mp, int bi, int bj) const {
    auto outerNSpace = nspace(cutoff, dim-1);
    auto ns = nspace(cutoff, 1);

    arma::uvec outerDims = arma::linspace<arma::uvec>(0, dim-1, dim);
    outerDims.shed_row(inner);

    arma::cx_double el = 0;

    #pragma omp parallel for reduction(compadd:el)
    for (std::size_t j = 0; j < outerNSpace.size(); ++j) {
      NPoint n(dim);
      n(outerDims) = outerNSpace[j];
      NPoint np = n;

      for (std::size_t ni = 0; ni < ns.size(); ++ni) {
        n(inner) = ns[ni](0);
        for (std::size_t nj = 0; nj < ns.size(); ++nj) {
          np(inner) = ns[nj](0);

          int fexp = (np(inner) - n(inner)) * N + (mp(inner) - m(inner));
          if (fexp != 0) {
            arma::cx_double phase = std::pow(-1, 1 + fexp);
            if (N%2 == 0) {
              phase *= std::exp(std::complex<double>(0, fexp*1.0/N*pi));
            }
            el += std::conj(_bs->blochC(m, n, bands[bi])) * _bs->blochC(mp, np, bands[bj]) *
              phase * std::complex<double>(0, 1.0/(2*pi*fexp));
          }
        }
      }
    }

    return el;
  }

}
