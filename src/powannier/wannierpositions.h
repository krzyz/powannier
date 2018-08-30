#ifndef POWANNIER_WANNIERPOSITIONS_H
#define POWANNIER_WANNIERPOSITIONS_H

#include <vector>

namespace POWannier {
  class WannierPositions {
    public:
      WannierPositions() {};

      WannierPositions(
          std::initializer_list<int> childrenNs);

      WannierPositions(
          std::vector<int> childrenNs);

      WannierPositions(
          std::initializer_list<WannierPositions> children);

      void addChild(WannierPositions child);

      int childrenNumber() const;
      int descendantsNumber() const;
      int descendantsNumberLeftTo(int i) const;

      WannierPositions getChild(int i);

    private:
      template <class It>
      WannierPositions(It first, It last);

      std::vector<WannierPositions> _children;
  };
}

#endif