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

  WannierPositions::WannierPositions(std::initializer_list<int> childrenNs) : 
      WannierPositions(childrenNs.begin(), childrenNs.end()) {}

  WannierPositions::WannierPositions(std::vector<int> childrenNs) : 
      WannierPositions(childrenNs.begin(), childrenNs.end()) {}


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
  
  int WannierPositions::descendantsNumber() const {
    if (_children.size() == 0) {
      return 1;
    }
    int value = 0;
    for (auto child : _children) {
      value += child.descendantsNumber();
    }
    return value;
  }

  int WannierPositions::descendantsNumberLeftTo(int i) const{
    if (_children.size() == 0) {
      return 0;
    }
    if (i >= _children.size()) {
      throw std::runtime_error("Attempting to reference nonexistant child!");
    }
    int value = 0;
    for (int j = 0; j < i; ++j) {
      value += _children[j].descendantsNumber();
    }
    return value;
  }

  WannierPositions WannierPositions::getChild(int i) {
    if (_children.size() == 0) {
      return *this;
    }
    if (i >= _children.size()) {
      throw std::runtime_error("Attempting to reference nonexistant child!");
    }
    return _children[i];
  }
}