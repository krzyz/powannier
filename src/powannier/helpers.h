#ifndef POWANNIER_HELPERS_H
#define POWANNIER_HELPERS_H

#include <complex>
#include <vector>

#include "aliases.h"

namespace POWannier {
  double chop(double number, double prec = 1e-15);
  std::complex<double> chop(std::complex<double> number, double prec = 1e-15);

  std::vector<NPoint> nspace(int cutoff, int dim);
}

#endif