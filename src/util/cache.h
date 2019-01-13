#ifndef _POWANNIER_CACHE_H
#define _POWANNIER_CACHE_H

#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

namespace POWannier {
  template <class V>
  class LazyInit {
    public:
      template<class Function>
      const V& get(Function&& initFunction) {
        std::call_once(_flag, [&] () {
            _data = std::make_unique<V>(initFunction());
          }   
        );  

        return *_data;
      }

    private:
      mutable std::once_flag _flag;
      mutable std::unique_ptr<V> _data;
  };

  template<class K, class V>
  class Cache {
    public:
      Cache() = default;

      template <class Function>
      Cache(Function&& initFunction) :
        _initFunction(initFunction) {}

      const V& get(K key) {
        {
          std::lock_guard<std::mutex> guard(_mutex);
          if (_cache.find(key) == _cache.end()) {
            _cache.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(key),
                std::forward_as_tuple()
            );  
          }   
        }

        return _cache[key].get(std::bind(_initFunction, key));
      }

    private:
      std::mutex _mutex;
      std::function<V(K)> _initFunction;
      std::map<K,LazyInit<V>> _cache;
  };
}

#endif