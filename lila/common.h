#pragma once

#include <complex>
#include <cstdint>

namespace lila {
using lila_size_t = int64_t;
using scomplex = std::complex<float>;
using complex = std::complex<double>;

template <class T> struct is_complex_t : public std::false_type {};
template <class T>
struct is_complex_t<std::complex<T>> : public std::true_type {};
template <class T> constexpr bool is_complex() {
  return is_complex_t<T>::value;
}
template <class T> constexpr bool is_real() { return !is_complex_t<T>::value; }

} // namespace lila
