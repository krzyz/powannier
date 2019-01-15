#include <functional>
#include <iostream>
#include <stdexcept>
#include "wannier.h"
#include "helpers.h"

namespace POWannier {
  Wannier::Wannier(std::shared_ptr<BlochSystem> bs,
                   const arma::cx_vec& coefficients,
                   std::vector<int> bands) :
    _bands(bands),
    _bs(std::move(bs)),
    _coefficients(coefficients) {

    Complex temp = 0;
    for (int i = 0; i < 10; ++i) {
      Position r;
      r.randn();
      temp += complexValue(r);
    }

    double beta = -atan(std::imag(temp)/real(temp));

    _coefficients *= std::exp(Complex(0,beta));
  }

  double Wannier::operator()(Position r) {
    return std::real(complexValue(r));
  }

  Complex Wannier::complexValue(Position r) {
    auto& ms = _bs->ms;

    arma::cx_double value = 0;
    #pragma omp parallel for reduction(compadd:value) collapse(2)
    for (std::size_t bandsi = 0; bandsi < _bands.size(); ++bandsi) {
      for (std::size_t mi = 0; mi < ms.size(); ++mi) {
        const auto& m = ms[mi];
        value += _coefficients[bandsi* ms.size() + mi] * _bs->bloch(m, r, _bands[bandsi]);
      }
    }

    return std::sqrt(1.0/(std::pow(_bs->N, _bs->dim()) * _bs->elementaryCellVolume())) * value;
  }

  void Wannier::draw(std::string fileName, int density, Position beg, Position end) {
    POWannier::draw(_bs->latticeBasis(), [&] (Position r) { return Position({operator()(r)});}, fileName, density,
      beg, end);
  }

  void Wannier::draw(std::string fileName, int density, DrawRange drawRange) {
    POWannier::draw(*_bs, [&] (Position r) { return Position({operator()(r)});}, fileName, density, drawRange);
  }

  const std::vector<int>& Wannier::bands() {
    return _bands;
  }
}
