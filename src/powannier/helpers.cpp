#include "helpers.h"

namespace POWannier {
  double chop(double number, double prec) {
    if (std::abs(number) < prec) {
      number = 0;
    }
    return number;
  }

  std::complex<double> chop(std::complex<double> number, double prec) {
    return std::complex<double>(chop(std::real(number), prec), chop(std::imag(number), prec));
  }

  std::vector<NPoint> space(int beg, int end, int dim) {
    std::vector<NPoint> ret;
    int cdim = end - beg + 1;
    int numberAll = pow(cdim, dim);
    for (int i = 0; i < numberAll; ++i) {
      int ih = i;
      NPoint n(dim);

      for (int j = dim-1; j >= 0; --j) {
        n(j) = ih % cdim + beg;
        ih /= cdim;
      }
      ret.push_back(n);
    }

    return ret;
  }

  std::vector<NPoint> nspace(int cutoff, int dim) {
    return space(-cutoff, cutoff, dim);
  }

  std::vector<NPoint> mspace(int N, int dim) {
    return space(0, N-1, dim);
  }

  int getIndex(const NPoint& n, int beg, int end) {
    int dim = n.n_elem;
    int ret = 0;
    int cdim = end - beg + 1;
    for (int j = 0; j < dim; ++j) {
      ret *= cdim;
      ret += n(j) - beg;
    }

    return ret;
  }

  int nIndex(const NPoint& n, int cutoff) {
    return getIndex(n, -cutoff, cutoff);
  }

  int mIndex(const NPoint& m, int N) {
    return getIndex(m, 0, N-1);
  }

  int bmIndex(const NPoint& m, int N, int bandi) {
    return mIndex(m, N) + bandi * std::pow(N, m.n_elem);
  }

}
