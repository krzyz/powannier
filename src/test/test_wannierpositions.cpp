#include <catch.hpp>
#include <iostream>

#include "../powannier/wannierpositions.h"

TEST_CASE("Check WannierPositions class", "[wannierpositions]") {
  SECTION("Set 3D Wannier Positions") {
    POWannier::WannierPositions pos({{1, 2}, {3}, {1}, {1, 1, 1, 1}});

    pos = pos.getChild(0);
  }
}
