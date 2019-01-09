#ifndef POWANNIER_BLOCHSYSTEMBUILDER_H
#define POWANNIER_BLOCHSYSTEMBUILDER_H

#include "potential.h"

namespace POWannier {
  class BlochSystem;

  class BlochSystemBuilder {
    public:
      BlochSystemBuilder(const Potential& potential) :
        _potential(potential),
        _laserWavelength(1),
        _fourierCutoff(2),
        _systemSize(5),
        _potentialDepth(8) {}

      BlochSystemBuilder& setLaserWavelength(double laserWavelength) {
        _laserWavelength = laserWavelength;
      }

      BlochSystemBuilder& setFourierCutoff(int fourierCutoff) {
        _fourierCutoff = fourierCutoff;
      }

      BlochSystemBuilder& setSystemSize(int systemSize) {
        _systemSize = systemSize;
      }

      BlochSystemBuilder& setPotentialDepth(int potentialDepth) {
        _potentialDepth = potentialDepth;
      }

      BlochSystem complete();

    private:
      const Potential _potential;
      double _laserWavelength;
      int _fourierCutoff;
      int _systemSize;
      int _potentialDepth;
  };
}

#endif