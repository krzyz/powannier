#ifndef POWANNIER_BLOCHSYSTEMBUILDER_H
#define POWANNIER_BLOCHSYSTEMBUILDER_H

#include "potential.h"

namespace POWannier {
  class BlochSystem;

  class BlochSystemBuilder {
    public:
      BlochSystemBuilder(const Potential& potential) :
        _potential(potential) {}

      BlochSystemBuilder& setLaserWavelength(double laserWavelength) {
        _laserWavelength = laserWavelength;

        return *this;
      }

      BlochSystemBuilder& setFourierCutoff(int fourierCutoff) {
        _fourierCutoff = fourierCutoff;

        return *this;
      }

      BlochSystemBuilder& setSystemSize(int systemSize) {
        _systemSize = systemSize;

        return *this;
      }

      BlochSystemBuilder& setPotentialDepth(int potentialDepth) {
        _potentialDepth = potentialDepth;

        return *this;
      }

      BlochSystem complete();

    private:
      const Potential _potential;
      double _laserWavelength = 1;
      int _fourierCutoff = 2;
      int _systemSize = 5;
      int _potentialDepth = 8;
  };
}

#endif
