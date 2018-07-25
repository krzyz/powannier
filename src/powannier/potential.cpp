#include "potential.h"
#include "../util/cubatureintegration.h"

namespace POWannier {
  using Potential = PotentialT<CubatureIntegration>;
}