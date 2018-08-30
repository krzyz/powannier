#include <stdexcept>
#include "rsystem.h"
#include "helpers.h"

namespace POWannier {
  RSystem::RSystem(std::shared_ptr<BlochSystem> bs, std::vector<int> bands) : 
    N(bs->N),
    dim(bs->dim),
    mdim(std::pow(N, dim)),
    bmdim(std::pow(N, dim) * bands.size()),
    cutoff(bs->cutoff),
    bands(bands),
    _bs(std::move(bs)) {
    std::vector<int> positions(bands.size(), 1);
    _wannierPositions = WannierPositions(positions);
    _r1EigenSystem = getREigensystem(0);
    std::cout << "r1 eigenvalues: " << std::endl << r1Eigenvalues() << std::endl;
  }

  RSystem::RSystem(std::shared_ptr<BlochSystem> bs, std::vector<int> bands, WannierPositions positions) : 
    N(bs->N),
    dim(bs->dim),
    mdim(std::pow(N, dim)),
    bmdim(std::pow(N, dim) * bands.size()),
    cutoff(bs->cutoff),
    bands(bands),
    _wannierPositions(positions),
    _bs(std::move(bs)) {
    if (_wannierPositions.descendantsNumber() != bands.size()) {
      throw std::runtime_error("Number of specified positions of wannier functions must be equal to bands number!");
    }
    _r1EigenSystem = getREigensystem(0);
    std::cout << "r1 eigenvalues: " << std::endl << r1Eigenvalues() << std::endl;
  }

  Wannier RSystem::getWannier(NPoint n) {
    NPoint elCellPositions(dim, arma::fill::zeros);
    return getWannier(n, elCellPositions);
  }


  Wannier RSystem::getWannier(NPoint n, NPoint elCellPositions) {
    if (n.n_elem != dim) {
      throw std::runtime_error("Specification of Wannier location must have the same dimension as RSystem!");
    }
    NPoint pos = n;
    if (N%2 == 0) {
      pos.for_each([&] (auto& x) {x += N/2 - 1;});
    } else {
      pos.for_each([&] (auto& x) {x += (N-1) / 2;});
    }

    arma::cx_mat eigenvectors = r1Eigenvectors();
    arma::vec eigenvalues = r1Eigenvalues();

    arma::cx_mat subspace(bmdim, bmdim, arma::fill::eye);
    arma::cx_mat subspace1(bmdim, bmdim, arma::fill::eye);

    auto currentWannierPositions = _wannierPositions;

    for (int i = 0; i < dim-1; ++i) {
      int idim = std::pow(N, dim-1-i);
      int currentPosInCell = elCellPositions(i);
      auto nextWannierPositions = currentWannierPositions.getChild(currentPosInCell);
      int subspaceStart = idim * bands.size() * pos(i) +
        currentWannierPositions.descendantsNumberLeftTo(currentPosInCell);
      int subspaceEnd = subspaceStart +
        idim * nextWannierPositions.descendantsNumber() - 1;

      {
        EigenSystem subsystemb = getSubREigensystem(i+1, subspace1);
        std::cout << "r" << i+2 << "eigenvalues of whole space: " << std::endl;
        std::cout << std::get<0>(subsystemb) << std::endl;
      }


      std::cout << "subspace eigenvalues: " << std::endl << eigenvalues.subvec(subspaceStart, subspaceEnd);
      subspace = subspace * eigenvectors.cols(subspaceStart, subspaceEnd);


      {
        EigenSystem subsystemb = getSubREigensystem(i, subspace);
        std::cout << "r" << i+1 << "eigenvalues inside subspace: " << std::endl;
        std::cout << std::get<0>(subsystemb) << std::endl;
      }

      EigenSystem subsystem = getSubREigensystem(i+1, subspace);

      eigenvalues = std::move(std::get<0>(subsystem));

      std::cout << "r" << i+2 << " eigenvalues: " << std::endl << eigenvalues;
      eigenvectors = std::move(std::get<1>(subsystem));
      currentWannierPositions = nextWannierPositions;
    }

    int wannierPosition = pos(dim-1) * bands.size() +
      currentWannierPositions.descendantsNumberLeftTo(elCellPositions(dim-1));
    arma::cx_mat coefficients = subspace * eigenvectors.col(wannierPosition);

    return Wannier(_bs, coefficients, bands);
  }

  arma::uvec RSystem::transformToBm() const {
    int idim = N;
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

  arma::uvec RSystem::transformFromInner(int inner) const {
    if (inner < 0 || inner > dim-1) {
      throw std::runtime_error("Inner dimension index no in range of dimension indices!");
    }
    if (inner == dim-1) {
      return arma::linspace<arma::uvec>(0, mdim-1, mdim);
    }

    arma::uvec indices(mdim);


    for (int i = 0; i < mdim; ++i) {
      auto m = _bs->ms[i];
      auto mrearranged = m;
      mrearranged.subvec(inner, dim-2) = m.subvec(inner+1, dim-1);
      mrearranged(dim-1) = m(inner);
      auto inew = mIndex(mrearranged, N);
      indices(inew) = i;
    }

    return indices;
  }

  arma::uvec RSystem::transformFromInnerToBm(int inner) const {
    arma::uvec indices = transformToBm();
    arma::uvec indicesFromInner = transformFromInner(inner);
    for (int i = 0; i < bands.size(); ++i) {
      arma::uvec indicesInOneBand = indices.subvec(i*mdim, (i+1) * mdim - 1);
      indices.subvec(i*mdim, (i+1) * mdim - 1) = indicesInOneBand(indicesFromInner);
    }

    return indices;
  }

  const arma::cx_mat& RSystem::r1Eigenvectors() {
    return std::get<1>(_r1EigenSystem);
  }

  const arma::vec& RSystem::r1Eigenvalues() {
    return std::get<0>(_r1EigenSystem);
  }

  void RSystem::sortEigensystem(arma::vec& eigenvalues, arma::cx_mat& eigenvectors) const {
    arma::uvec sortedIndices = arma::sort_index(eigenvalues);
    eigenvalues = eigenvalues(sortedIndices);
    eigenvectors = eigenvectors.cols(sortedIndices);
  }

  EigenSystem RSystem::getREigensystem(int inner) const {
    int submdim = N*bands.size();

    auto outerMSpace = mspace(N, dim-1);
    auto omdim = outerMSpace.size();

    arma::uvec outerDims = arma::linspace<arma::uvec>(0, dim-1, dim);
    outerDims.shed_row(inner);

    arma::vec reigval(bmdim);
    arma::cx_mat reigvec(bmdim, bmdim, arma::fill::zeros);

    for (int i = 0; i < omdim; ++i) {

      NPoint m(dim);
      m(outerDims) = outerMSpace[i];
      NPoint mp = m;

      arma::cx_mat r1dBands(submdim, submdim);

      for (int bi = 0; bi < bands.size(); ++bi) {
        int first_row = bi*N;
        for (int bj = bi; bj < bands.size(); ++bj) {
          int first_col = bj*N;
          std::cout << first_row << " " << first_col << std::endl;
          std::cout << arma::abs(rInnerMatrix(inner, m, bi, bj)) << std::endl;
          r1dBands.submat(first_row, first_col, arma::size(N, N)) = 
            rInnerMatrix(inner, m, bi, bj);
        }
      }

      arma::vec eigval;
      arma::cx_mat eigvec;

      std::cout << "r1dBands: " << std::endl;
      std::cout << arma::abs(r1dBands) << std::endl;

      arma::eig_sym(eigval, eigvec, r1dBands);

      reigval.subvec(i*submdim, arma::size(eigval)) = std::move(eigval);
      reigvec.submat(i*submdim, i*submdim, arma::size(eigvec)) = std::move(eigvec);
    }

    auto indicesInStandardBasis = transformFromInnerToBm(inner);
    std::cout << indicesInStandardBasis << std::endl;

    //reigval = reigval(indicesInStandardBasis);
    reigvec = reigvec.rows(indicesInStandardBasis);

    sortEigensystem(reigval, reigvec);

    return make_tuple(std::move(reigval), std::move(reigvec));
  }

  EigenSystem RSystem::getSubREigensystem(int inner, const arma::cx_mat& subspace) const {
    auto rEigensystem = getREigensystem(inner);
    const auto& rEigval = std::get<0>(rEigensystem);
    const auto& rEigvec = std::get<1>(rEigensystem);
    arma::cx_mat subr = 
      subspace.t() * 
      rEigvec * arma::diagmat(rEigval) * rEigvec.t() *
      subspace;

    std::cout << subr << std::endl;

    arma::vec eigval;
    arma::cx_mat eigvec;
    arma::eig_sym(eigval, eigvec, subr);

    return std::make_tuple(std::move(eigval), std::move(eigvec));
  }

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

  arma::cx_double RSystem::sumOverNs(int inner, NPoint m, NPoint mp, int bi, int bj) const {
    auto outerNSpace = nspace(cutoff, dim-1);
    auto ns = nspace(cutoff, 1);

    arma::uvec outerDims = arma::linspace<arma::uvec>(0, dim-1, dim);
    outerDims.shed_row(inner);

    arma::cx_double el = 0;

    for (int j = 0; j < outerNSpace.size(); ++j) {
      NPoint n(dim);
      n(outerDims) = outerNSpace[j];
      NPoint np = n;

      for (int ni = 0; ni < ns.size(); ++ni) {
        n(inner) = ns[ni](0);
        for (int nj = 0; nj < ns.size(); ++nj) {
          np(inner) = ns[nj](0);

          double fexp = (np(inner) - n(inner)) * N + (mp(inner) - m(inner));
          if (fexp != 0) {
            arma::cx_double phase = std::pow(-1, 1 + fexp);
            if (N%2 == 0) {
              phase *= std::exp(std::complex<double>(0, fexp/N*pi));
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
