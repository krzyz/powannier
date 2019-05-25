#ifndef POWANNIER_WANNIERPRODUCT_H
#define POWANNIER_WANNIERPRODUCT_H

#include "wannier.h"
#include "aliases.h"

namespace POWannier {
  class WannierProduct : public Wannier {
    public:
      static std::unique_ptr<WannierProduct> createUnique(std::initializer_list<std::shared_ptr<Wannier>> wanniers);

      double operator()(Position r) const override;

      arma::cx_double blochCoefficient(int mi, int bandi, int ni) const override;

      const std::vector<int>& bands() const override;

      const Wannier& get1dWannier(int n) const override;

      arma::vec energies(NPoint m) const override;

      std::vector<POWannier::NPoint> ms() const override;

      int cutoff() const override;

      int N() const override;

      int dim() const override {
        return _dim;
      }

    private:
      WannierProduct(std::shared_ptr<Wannier>> wanniers);
      std::vector<std::shared_ptr<Wannier>> _wanniers;
      int _dim;
  };
}

#endif