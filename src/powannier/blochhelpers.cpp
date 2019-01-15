#include "blochhelpers.h"
#include <iostream>

namespace POWannier {
  void drawBands(int numberOfBands, const Potential& V, double kl, int cutoff, double s, std::string fileName, int density, DrawRange drawRange) {
    auto bs = BlochSystem::createWithPotential(V)
      .setLaserWavelength(kl)
      .setFourierCutoff(cutoff)
      .setPotentialDepth(s)
      .complete();

    ReciprocalVector beg, end;
    if (drawRange == DrawRange::BrillouinZone) {
      beg = -pi * ReciprocalVector(bs.dim(), arma::fill::ones);
      end = pi * ReciprocalVector(bs.dim(), arma::fill::ones);

      draw(bs.latticeBasis(), [&] (ReciprocalVector k) -> arma::vec {
        auto bloch = Bloch(k, V, kl, cutoff, s);
        return bloch.energies().head(numberOfBands);
      }, fileName, density, beg, end);
    } else {
      throw std::runtime_error("invalid draw range!");
    }
  }
}
