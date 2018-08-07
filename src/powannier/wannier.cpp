#include "wannier.h"

namespace POWannier {
  Wannier::Wannier(std::shared_ptr<BlochSystem> bs,
                   const arma::cx_vec& coefficients) :
    _bs(std::move(bs)),
    _coefficients(coefficients) {}

  double Wannier::operator()(Position r) {
    arma::cx_double value;

    auto& ms = _bs->ms;

    for (size_t mi = 0; mi < ms.size(); ++mi) {
      const auto& m = ms[mi];
      value += _coefficients[mi] * _bs->bloch(m, r);
    }

    return std::real(value);
  }

}