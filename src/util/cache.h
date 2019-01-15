#ifndef _POWANNIER_CACHE_H
#define _POWANNIER_CACHE_H

#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

namespace POWannier {

  /**
   * @brief Thread-safe cache that can be used to store 
   * the results of a given single argument function.
   * 
   * @details
   * Usage:
   * @code
   * Cache<int,double> cache([](int x) { return x/2.0; });
   * double y = cache.get(1); // y = 0.5
   * @endcode
   * @tparam K argument of the function.
   * @tparam V return type of the function.
   */
  template<class K, class V>
  class Cache {
    private:
      class ValueInit {
        public:
          template<class Function>
          const V& get(Function&& initFunction) {
            std::call_once(_flag, [&] () {
                _data = initFunction();
              }   
            );  

            return _data;
          }

        private:
          mutable std::once_flag _flag;
          mutable V _data;
      };

    public:
      Cache() = default;

      /**
       * Constructor.
       * 
       * @param initFunction A function whose return values are stored in the cache.
       * 
       * @pre @p initFunction must take a value of type @p K as its only
       * argument and return a value of type @p V. Currently @p K must 
       * support comparison (as an underlying data structure is @p std::map).
       */
      template <class Function>
      Cache(Function&& initFunction) :
        _initFunction(initFunction) {}

      /**
       * Get a cached value returned by the @p initFunction function for a given argument
       * (generate it first if it is the first time this value is needed).
       * 
       * @param key the argument for the function
       */
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
      std::map<K,ValueInit> _cache;
  };
}

#endif