#pragma once

#include "../common.h"
#include <complex>

#ifdef LILA_USE_MKL
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
#endif
#ifdef LILA_USE_ACCELERATE
#include <Accelerate/Accelerate.h>
#endif

namespace lila {
namespace blaslapack {

#ifdef LILA_USE_MKL // Using the Intel MKL
using blas_size_t = int;
using blas_float_t = float;
using blas_double_t = double;
using blas_scomplex_t = MKL_Complex8;
using blas_complex_t = MKL_Complex16;
using lapack_ret_t = void;
#define LAPACK_CONST const
#define __LAPACK_ROUTINE_NAME(x) x
#endif

#ifdef LILA_USE_LAPACK // Using normal LAPACK
using blas_size_t = int;
using blas_float_t = float;
using blas_double_t = double;
using blas_scomplex_t = std::complex<float>;
using blas_complex_t = std::complex<double>;
using lapack_ret_t = void;
#define LAPACK_CONST const
#define __LAPACK_ROUTINE_NAME(x) x##_
#endif

#ifdef LILA_USE_ACCELERATE // Using OSX Accellerate
using blas_size_t = int;
using blas_float_t = float;
using blas_double_t = double;
using blas_scomplex_t = __CLPK_complex;
using blas_complex_t = __CLPK_doublecomplex;
using lapack_ret_t = int;
#define LAPACK_CONST
#define __LAPACK_ROUTINE_NAME(x) x##_
#endif

// Define template alias to convert lila
// types to blas types
template <class T> struct blas_ts;
template <> struct blas_ts<int> { using type = blas_size_t; };
template <> struct blas_ts<float> { using type = blas_float_t; };
template <> struct blas_ts<double> { using type = blas_double_t; };
template <> struct blas_ts<std::complex<float>> {
  using type = blas_scomplex_t;
};
template <> struct blas_ts<std::complex<double>> {
  using type = blas_complex_t;
};

template <class T> using blas_t = typename blas_ts<T>::type;

// Define template alias to convert blas
// types to lila types
template <class T> struct lila_ts;
template <> struct lila_ts<blas_size_t> { using type = int; };
template <> struct lila_ts<blas_float_t> { using type = float; };
template <> struct lila_ts<blas_double_t> { using type = double; };
template <> struct lila_ts<blas_scomplex_t> {
  using type = std::complex<float>;
};
template <> struct lila_ts<blas_complex_t> {
  using type = std::complex<double>;
};

template <class T> using lila_t = typename lila_ts<T>::type;

#define LILA_BLAS_CAST(COEFFT, PTR)                                            \
  reinterpret_cast<blaslapack::blas_t<COEFFT> *>(PTR)
#define LILA_BLAS_CONST_CAST(COEFFT, PTR)                                      \
  reinterpret_cast<blaslapack::blas_t<COEFFT> *>(const_cast<COEFFT *>(PTR))

///////////////////////////////////////////
// Convert a blas number to a lila number
template <class T> inline lila_t<T> blas_to_lila(T const &);
template <> inline lila_t<blas_float_t> blas_to_lila(blas_float_t const &x) {
  return x;
}
template <> inline lila_t<blas_double_t> blas_to_lila(blas_double_t const &x) {
  return x;
}

#ifdef LILA_USE_MKL // Using the Intel MKL
template <>
inline lila_t<blas_scomplex_t> blas_to_lila(blas_scomplex_t const &x) {
  return x;
}
template <>
inline lila_t<blas_complex_t> blas_to_lila(blas_complex_t const &x) {
  return x;
}
#endif

#ifdef LILA_USE_LAPACK // Using normal LAPACK
template <>
inline lila_t<blas_scomplex_t> blas_to_lila(blas_scomplex_t const &x) {
  return x;
}
template <>
inline lila_t<blas_complex_t> blas_to_lila(blas_complex_t const &x) {
  return x;
}
#endif

#ifdef LILA_USE_ACCELERATE // Using OSX Accellerate
template <>
inline lila_t<blas_scomplex_t> blas_to_lila(blas_scomplex_t const &x) {
  return std::complex<float>(x.r, x.i);
}
template <>
inline lila_t<blas_complex_t> blas_to_lila(blas_complex_t const &x) {
  return std::complex<double>(x.r, x.i);
}
#endif

} // namespace blaslapack
} // namespace lila
