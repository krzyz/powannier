#include <functional>
#include <iostream>
#include <stdexcept>
#include "wannier.h"

namespace POWannier {
  Wannier::Wannier(std::shared_ptr<BlochSystem> bs,
                   const arma::cx_vec& coefficients,
                   std::vector<int> bands) :
    bands(bands),
    _bs(std::move(bs)),
    _coefficients(coefficients) {}

  double Wannier::operator()(Position r) {
    arma::cx_double value;

    auto& ms = _bs->ms;

    #pragma omp parallel for reduction(+:value) collapse(2)
    for (int bandsi = 0; bandsi < bands.size(); ++bandsi) {
      for (size_t mi = 0; mi < ms.size(); ++mi) {
        const auto& m = ms[mi];
        value += _coefficients[bandsi* ms.size() + mi] * _bs->bloch(m, r, bands[bandsi]);
      }
    }

    return std::sqrt(1.0/(std::pow(_bs->N, _bs->dim) * _bs->V->elementaryCellVolume())) * std::real(value);
  }

  void Wannier::drawB(std::string filename, int density, Position beg, Position end) {
    std::cout << beg << end << std::endl;
    //auto points = mspace(density, dim);
    arma::vec range = arma::linspace<arma::vec>(0, 1, density);
    std::ofstream outfile(filename);
    outfile.precision(14);

    std::function<void(Position,int)> nDloop;
    nDloop = [&] (Position alpha, int n) {
      if (n == _bs->dim) {
        Position r = alpha * _bs->V->latticeBasis();
        arma::rowvec toPrint(_bs->dim + 1);
        toPrint.head(_bs->dim) = r;
        toPrint.tail(1) = operator()(r);
        toPrint.raw_print(outfile);
      } else {
        for (const auto& x : range) {
          alpha(n) = beg(n) + x * (end(n) - beg(n));
          nDloop(alpha, n+1);
        }
      }
    };

    Position alpha(_bs->dim);
    nDloop(alpha, 0);
  }

  void Wannier::drawB(std::string filename, int density, DrawRange drawRange) {
    Position beg, end;
    if (drawRange == WholeLattice) {
      double ibeg = -_bs->N/2.0;
      double iend =  _bs->N/2.0;
      if (_bs->N %2 == 0) {
        ibeg = -(_bs->N - 1) / 2.0;
        iend =  (_bs->N + 1) / 2.0;
      }
      beg = ibeg * Position(_bs->dim, arma::fill::ones);
      end = iend * Position(_bs->dim, arma::fill::ones);
      drawB(filename, density, beg, end);
    } else if (drawRange == ElementaryCell) {
      beg = Position(_bs->dim, arma::fill::zeros);
      end = Position(_bs->dim, arma::fill::ones);
      drawB(filename, density, beg, end);
    } else {
      throw std::runtime_error("invalid draw range!");
    }
  }
}
