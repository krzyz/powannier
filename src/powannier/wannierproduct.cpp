#include "wannierproduct.h"

namespace POWannier {
  WannierProduct::WannierProduct(std::initializer_list<std::shared_ptr<Wannier>> wanniers) :
    _wanniers(wanniers),
    _dim(std::accumulate(std::begin(wanniers), std::end(wanniers), 0, [](int current, auto wannier) {
      return current + wannier->dim();
    })) {}

  std::unique_ptr<WannierProduct> createUnique(std::initializer_list<std::shared_ptr<Wannier>> wanniers) {
    return std::make_unique(new WannierProduct(wanniers));
  }

  double WannierProduct::operator()(Position r) const override {
    double ret = 1;

    int dimInd = 0;

    for (auto wannier : _wanniers) {
      nextInd = dimInd + wannier->dim();

      ret *= (*wannier)(r.subvec(dimInd, nextInd-1);
      dimInd = nextInd;
    }

    return ret;
  }

  arma::cx_double WannierProduct::blochCoefficient(int mi, int bandi, int ni) const override;

  const std::vector<int>& WannierProduct::bands() const override;

  const Wannier& get1dWannier(int n) const override;
}