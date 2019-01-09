#include "blochsystem.h"
#include "blochsystembuilder.h"

namespace POWannier {
  BlochSystem BlochSystemBuilder::complete() {
    return BlochSystem(_potential, _laserWavelength, _fourierCutoff, _systemSize, _potentialDepth);
  }
}