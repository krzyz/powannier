#include <functional>
#include <iostream>
#include <stdexcept>
#include "wannier.h"
#include "helpers.h"

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

    return std::sqrt(1.0/(std::pow(_bs->N, _bs->dim()) * _bs->elementaryCellVolume())) * std::real(value);
  }

  void Wannier::drawB(std::string fileName, int density, Position beg, Position end) {
    draw(*_bs, [&] (Position r) { return Position({operator()(r)});}, fileName, density,
      beg, end);
  }

  void Wannier::drawB(std::string fileName, int density, DrawRange drawRange) {
    draw(*_bs, [&] (Position r) { return Position({operator()(r)});}, fileName, density, drawRange);
  }
}
