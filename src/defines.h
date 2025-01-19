//
// Created by Andrew Astakhov on 18.01.25.
//

#ifndef DEFINES_H
#define DEFINES_H

#include <limits>
#include <cmath>

template<typename T>
constexpr static T CMMCORE_EPSILON() noexcept {
  return std::numeric_limits<T>::epsilon();
}
template<typename T>
constexpr static T CMMCORE_SQRT_EPSILON() noexcept {
  return std::sqrt(std::numeric_limits<T>::epsilon());
}
template<typename T>
constexpr static T CMMCORE_MAX() noexcept {
  return std::sqrt(std::numeric_limits<T>::max());
}
template<typename T>
constexpr static T CMMCORE_MIN() noexcept {
  return std::sqrt(std::numeric_limits<T>::min());
}
#endif //DEFINES_H
