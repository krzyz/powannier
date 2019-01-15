#include <stdexcept>
#include "wannierpositions.h"

namespace POWannier {
  template <class It>
  WannierPositions::WannierPositions(It first, It last) {
    for (auto it = first; it != last; ++it) {
      WannierPositions child;
      auto n = *it;
      for (int i = 0; i < n; ++i) {
        child.addChild(WannierPositions());
      }
      addChild(child);
    }
  }

  WannierPositions::WannierPositions(std::initializer_list<int> leafNs) : 
      WannierPositions(leafNs.begin(), leafNs.end()) {}

  WannierPositions::WannierPositions(std::vector<int> leafNs) : 
      WannierPositions(leafNs.begin(), leafNs.end()) {}


  WannierPositions::WannierPositions(
      std::initializer_list<WannierPositions> children) {
    for (auto child : children) {
      addChild(child);
    }
  }

  void WannierPositions::addChild(WannierPositions child) {
    _children.push_back(child);
  }

  int WannierPositions::childrenNumber() const {
    return _children.size();
  }
  
  int WannierPositions::descendantLeavesNumber() const {
    if (_children.size() == 0) {
      return 1;
    }
    int value = 0;
    for (auto child : _children) {
      value += child.descendantLeavesNumber();
    }
    return value;
  }

  int WannierPositions::descendantLeavesNumberLeftTo(int i) const{
    if (_children.size() == 0) {
      return 0;
    }
    if (i >= static_cast<int>(_children.size())) {
      throw std::runtime_error("Attempting to reference nonexistant child!");
    }
    int value = 0;
    for (int j = 0; j < i; ++j) {
      value += _children[j].descendantLeavesNumber();
    }
    return value;
  }

  WannierPositions WannierPositions::getChild(int i) {
    if (_children.size() == 0) {
      return *this;
    }
    if (i >= static_cast<int>(_children.size())) {
      throw std::runtime_error("Attempting to reference nonexistant child!");
    }
    return _children[i];
  }
}
