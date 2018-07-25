#include "helpers.h"

namespace POWannier {
  double chop(double number, double prec) {
    if (abs(number) < prec)
      number = 0;
    return number;
  }

  std::complex<double> chop(std::complex<double> number, double prec) {
    return std::complex<double>(chop(std::real(number)), chop(std::imag(number)));
  }

  std::vector<NPoint> nspace(int cutoff, int dim) {
    std::vector<NPoint> ret;
    int number1D = 2 * cutoff + 1;
    int numberAll = pow(number1D, dim);
    for (int i = 0; i < numberAll; ++i) {
      int ih = i;
      NPoint n(dim);
      for (int j = 0; j < dim; ++j) {
        n(j) = ih % number1D - cutoff;
        ih /= number1D;
      }
      ret.push_back(n);
    }

    return ret;
  }
}