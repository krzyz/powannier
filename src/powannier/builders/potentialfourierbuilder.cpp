#include "potentialfourierbuilder.h"
#include "potentialfourierbuilderaddition.h"

namespace POWannier {
  PotentialFourierBuilderAddition PotentialFourierBuilder::provideMultiple(std::vector<NPoint> indices) {
    std::queue<NPoint> indicesQ;

    for (auto index : indices) {
      indicesQ.push(index);
    }

    return PotentialFourierBuilderAddition{_potential, indicesQ};
  }

  PotentialFourierBuilderAddition PotentialFourierBuilder::provideOne(NPoint index) {
    std::queue<NPoint> indices;
    indices.push(index);
    return PotentialFourierBuilderAddition{_potential, indices};
  }
}