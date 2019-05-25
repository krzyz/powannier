#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
#include "wannier.h"
#include "helpers.h"

namespace POWannier {
  class InteractionContext {
    using WannierCRef = std::reference_wrapper<const Wannier>;
    using TwoWanniersWithQ = std::tuple<WannierCRef,WannierCRef,int,int,int>;
    using FourWanniers = std::tuple<WannierCRef,WannierCRef,WannierCRef,WannierCRef,int,int,int,int>;

    static Cache<TwoWanniersWithQ,arma::cx_double> integration2Cache;
    static Cache<FourWanniers,arma::cx_double> integration4Cache;

    arma::cx_double calculate1dIntegral2(TwoWanniersWithQ twoWanniersWithQ) {
      auto& wan1 = std::get<0>(twoWanniersWithQ);
      auto& wan2 = std::get<1>(twoWanniersWithQ);
      int i = std::get<2>(twoWanniersWithQ);
      int j = std::get<3>(twoWanniersWithQ);
      int qi = std::get<4>(twoWanniersWithQ);
      
      auto& ms1 = wan1.ms();
      auto& ms2 = wan2.ms();

      int N = wan1.N();

      // ensureSameSize(wan1, wan2);

      arma::cx_double value = 0;

      #pragma omp parallel for reduction(compadd:value) collapse(2)
      for (std::size_t m1i = 0; m1i < ms1.size(); ++m1i) {
        for (std::size_t band1i = 0; band1i < wan1.bands().size(); ++band1i) {
          arma::cx_double part1 = 0;

          for (std::size_t m2i = 0; m2i < ms2.size(); ++m2i) {
            if ((m1i - m2i + qi) % N == 0) {
              for (std::size_t band2i = 0; band2i < wan2.bands().size(); ++band2i) {
                arma::cx_double part2 = 0;

                for (auto n1 : wan1.BlochNs()) {
                  int n2 = (m1i - m2i + qi) / N + n1;

                  part2 += conj(wan1.blochCoefficient(m1i, band2i, n1)) *
                                wan2.blochCoefficient(m2i, band2i, n2);
                }
              }
              part1 += part2 * wan2.coefficient(m2i, band2i) * std::exp(std::complex<double>(0, -2*pi*j*m2i/(double)N));
            }
          }
          value += part1 * conj(wan2.coefficient(m1i, band1i)) * std::exp(std::complex<double>(0, 2*pi*i*m1i/(double)N));
        }
      }

      return value;
    }

    template <class Function>
    static arma::cx_double calculateOuterIntegral(const Wannier& wan1di, const Wannier& wan1dj, const Wannier& wan1dk, const Wannier& wan1dl,
      int i, int j, int k, int l, Function&& innerFunction) {

      complex<double> value = 0;

      int qs = 2 * wan1di.N() * (wan1di.cutoff() + 0.5) - 1;

      for (int q1d = -qs; q1d <= qs; ++q1d) {

        TwoWanniersWithQ wanil = std::make_tuple(std::cref(wan1di), std::cref(wan1dl), i, l, -q1d);
        TwoWanniersWithQ wanjk = std::make_tuple(std::cref(wan1dj), std::cref(wan1dk), j, k, q1d);

        value += innerFunction(q1d) * 
          integration2Cache.get(wanil) * integration2Cache.get(wanjk);
      }

      return value;
    }

    static arma::cx_double calculate1dIntegral4(FourWanniers fourWanniers) {
      Wannier& wan1di = std::get<0>(fourWanniers);
      Wannier& wan1dj = std::get<1>(fourWanniers);
      Wannier& wan1dk = std::get<2>(fourWanniers);
      Wannier& wan1dl = std::get<3>(fourWanniers);
      int i = std::get<4>(fourWanniers);
      int j = std::get<5>(fourWanniers);
      int k = std::get<6>(fourWanniers);
      int l = std::get<7>(fourWanniers);

      auto innerFunction = [] (int) {
        return 1;
      };

      return calculateOuterIntegral(wan1di, wan1dj, wan1dk, wan1dl, i, j, k, l, innerFunction);
    }

    template <class Function>
    static arma::cx_double get1dIntegral(const Wannier& wani, const Wannier& wanj, const Wannier& wank, const Wannier& wanl,
      NPoint ip, NPoint jp, NPoint kp, NPoint lp, Function&& interactionFunction, std::vector<int> funDepDims, NPoint q) {
      int n = q._dim;

      const Wannier& wan1di = wani.getWannier1d(n);
      const Wannier& wan1dj = wanj.getWannier1d(n);
      const Wannier& wan1dk = wank.getWannier1d(n);
      const Wannier& wan1dl = wanl.getWannier1d(n);
      int i = ip(n);
      int j = jp(n);
      int k = kp(n);
      int l = lp(n);

      qh.head(n) = q;

      if (std::find(funDepDims.begin(), funDepDims.end(), n) == funDepDims.end()) {
        qh.tail(1) = 0;

        FourWanniers wanijkl = std::make_tuple(std::cref(wan1di), std::cref(wan1dj),
          std::cref(wan1dk), std::cref(wan1dl), i, j, k, l);

        return integration4Cache.get(wanijkl) * get1dIntegral(wani, wanj, wank, wanl, ip, jp, kp, lp, qh);
      } else {
        auto innerFunction = [] (int q1d) {
          qh.tail(1) = q1d;

          if (n == wani.dim()-1) {
            return interactionFunction(2 * pi * qh / wani.N());
          } else {
            return get1dIntegral(wani, wanj, wank, wanl, i, j, k, l, interactionFunction, qh);
          }
        };

        return calculateOuterIntegral(wan1di, wan1dj, wan1dk, wan1dl, i, j, k, l, innerFunction);
      }
    }
  };

  static double tunneling(const Wannier& wan1, const Wannier& wan2, NPoint l) {
    // ensureCompatible(wan1, wan2);

    auto& ms = wan1.ms();

    arma::cx_double value = 0;

    #pragma omp parallel for reduction(compadd:value) collapse(2)
    for (auto band : commonBands(wan1, wan2)) {
      auto band1i = wan1.getLocalBandIndex(band);
      auto band2i = wan2.getLocalBandIndex(band);

      for (std::size_t mi = 0; mi < ms.size(); ++mi) {
        const auto& m = ms[mi];
        arma::cx_double el = conj(wan1.coefficient(mi, band1i)) * wan2.coefficient(mi, band2i);
        el *= wan1.energies(m)[band] * std::exp(std::complex<double>(0, -2 * arma::dot(l, m) * pi / wan1.N()));
        value += el;
      }
    }

    return value;
  }

  template <class Function>
  static arma::cx_double interaction(
    const Wannier& wani, const Wannier& wanj, const Wannier& wank, const Wannier& wanl,
    NPoint ip, NPoint jp, NPoint kp, NPoint lp) {
    
    std::vector<int> funDepDims(0);
    auto interactionFunction = [] (NPoint k) {
      return 1;
    }

    return interaction(wani, wanj, wank, wanl, ip, jp, kp, lp, interactionFunction, funDepDims, q)
  }

  template <class Function>
  static arma::cx_double interaction(
    const Wannier& wani, const Wannier& wanj, const Wannier& wank, const Wannier& wanl,
    NPoint ip, NPoint jp, NPoint kp, NPoint lp, Function&& interactionFunction) {

    std::vector<int> funDepDims(wani.dim());
    std::iota(funDepDims.begin(), funDepDims.end(), 0);

    return interaction(wani, wanj, wank, wanl, ip, jp, kp, lp, interactionFunction, funDepDims, q)
  }

  template <class Function>
  static arma::cx_double interaction(
    const Wannier& wani, const Wannier& wanj, const Wannier& wank, const Wannier& wanl,
    NPoint ip, NPoint jp, NPoint kp, NPoint lp, Function&& interactionFunction, std::vector<int> funDepDims) {

    using InteractionContext;
    Cache<WanWithInd,arma::cx_double> integrationCache(calculate1dIntegral2);
    Cache<K,arma::cx_double> integration4Cache(calculate1dIntegral4);

    NPoint q(0);

    return get1dIntegral(wani, wanj, wank, wanl, ip, jp, kp, lp, interactionFunction, funDepDims, q)
  }


}