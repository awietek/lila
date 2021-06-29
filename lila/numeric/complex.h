#pragma once

#include <complex>

#include <lila/common.h>
#include <lila/detail/complex_detail.h>
#include <lila/matrix.h>
#include <lila/vector.h>

namespace lila {

template <class coeff_t>
using real_t = typename detail::real_type_struct<coeff_t>::type;

template <class coeff_t>
using complex_t = typename detail::complex_type_struct<coeff_t>::type;

inline float real(float x) { return x; }
inline double real(double x) { return x; }
inline float real(std::complex<float> x) { return x.real(); }
inline double real(std::complex<double> x) { return x.real(); }

inline float imag(float) { return 0.; }
inline double imag(double) { return 0.; }
inline float imag(std::complex<float> x) { return x.imag(); }
inline double imag(std::complex<double> x) { return x.imag(); }

inline float conj(float x) { return x; }
inline double conj(double x) { return x; }
inline std::complex<float> conj(std::complex<float> x) { return std::conj(x); }
inline std::complex<double> conj(std::complex<double> x) {
  return std::conj(x);
}

template <class coeff_t>
inline Matrix<real_t<coeff_t>> Real(Matrix<coeff_t> const &X) {
  Matrix<real_t<coeff_t>> Y(X.nrows(), X.ncols());
  for (int i = 0; i < X.nrows(); ++i)
    for (int j = 0; j < X.ncols(); ++j)
      Y(i, j) = lila::real(X(i, j));
  return Y;
}

template <class coeff_t>
inline Vector<real_t<coeff_t>> Real(Vector<coeff_t> const &X) {
  Vector<real_t<coeff_t>> Y(X.nrows());
  for (int i = 0; i < X.nrows(); ++i)
    Y(i) = lila::real(X(i));
  return Y;
}

template <class coeff_t>
inline Matrix<real_t<coeff_t>> Imag(Matrix<coeff_t> const &X) {
  Matrix<real_t<coeff_t>> Y(X.nrows(), X.ncols());
  for (int i = 0; i < X.nrows(); ++i)
    for (int j = 0; j < X.ncols(); ++j)
      Y(i, j) = lila::imag(X(i, j));
  return Y;
}

template <class coeff_t>
inline Vector<real_t<coeff_t>> Imag(Vector<coeff_t> const &X) {
  Vector<real_t<coeff_t>> Y(X.nrows());
  for (int i = 0; i < X.nrows(); ++i)
    Y(i) = lila::imag(X(i));
  return Y;
}

template <class coeff_t> Matrix<coeff_t> Conj(Matrix<coeff_t> const &X) {
  Matrix<coeff_t> X_c(X.nrows(), X.ncols());
  for (int i = 0; i < X.nrows(); ++i)
    for (int j = 0; j < X.ncols(); ++j)
      X_c(i, j) = lila::conj(X(i, j));
  return X_c;
}

template <class coeff_t> Matrix<coeff_t> Herm(Matrix<coeff_t> const &X) {
  return Conj(Transpose(X));
}

} // namespace lila
