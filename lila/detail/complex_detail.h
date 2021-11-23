#pragma once
#include <complex>

namespace lila::detail {

template <class coeff_t> struct real_type_struct { typedef coeff_t type; };

template <class coeff_t> struct real_type_struct<std::complex<coeff_t>> {
  typedef coeff_t type;
};

template <class coeff_t> struct complex_type_struct {
  typedef std::complex<coeff_t> type;
};

template <class coeff_t> struct complex_type_struct<std::complex<coeff_t>> {
  typedef std::complex<coeff_t> type;
};

} // namespace lila::detail
