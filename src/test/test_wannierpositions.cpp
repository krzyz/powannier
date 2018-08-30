#include <catch.hpp>
#include <iostream>

#include "../powannier/wannierpositions.h"

TEST_CASE("Check WannierPositions class", "[wannierpositions]") {
  SECTION("Set 3D Wannier Positions") {
    POWannier::WannierPositions pos({{1, 2}, {3}, {1}, {1, 1, 1, 1}});

    std::cout << pos.descendantsNumber() << std::endl;
    std::cout << pos.descendantsNumberLeftTo(2) << std::endl;
    pos = pos.getChild(0);
    std::cout << pos.descendantsNumberLeftTo(1) << std::endl;
  }
}
