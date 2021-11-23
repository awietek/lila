#pragma once

#include <cassert>
#include <complex>

#include <lila/matrix.h>
#include <lila/vector.h>

#include <lila/blaslapack/blaslapack.h>

namespace lila {

template <class coeff_t> inline void Scale(coeff_t alpha, Vector<coeff_t> &v) {
  blas_size_t size = v.size();
  blas_size_t inc = 1;
  blaslapack::scal(&size, LILA_BLAS_CONST_CAST(coeff_t, &alpha),
                   LILA_BLAS_CAST(coeff_t, v.data()), &inc);
}

template <class coeff_t>
inline void Scale(coeff_t alpha, VectorView<coeff_t> v) {
  blas_size_t size = v.size();
  blas_size_t inc = v.inc();
  blaslapack::scal(&size, LILA_BLAS_CONST_CAST(coeff_t, &alpha),
                   LILA_BLAS_CAST(coeff_t, v.data()), &inc);
}

template <class coeff_t> inline void Scale(coeff_t alpha, Matrix<coeff_t> &A) {
  blas_size_t dx = A.size();
  blas_size_t inc = 1;
  blaslapack::scal(&dx, LILA_BLAS_CONST_CAST(coeff_t, &alpha),
                   LILA_BLAS_CAST(coeff_t, A.data()), &inc);
}

template <class coeff_t>
inline void Scale(coeff_t alpha, MatrixView<coeff_t> A) {
  blas_size_t ld = A.ld();
  blas_size_t m = A.m();
  blas_size_t n = A.n();
  blas_size_t incm = A.incm();
  blas_size_t incn = A.incn();

  // Perform a column-wise scale
  for (blas_size_t col = 0; col < n; ++col) {
    blaslapack::scal(&m, LILA_BLAS_CONST_CAST(coeff_t, &alpha),
                     LILA_BLAS_CAST(coeff_t, A.data() + col * incn * ld),
                     &incm);
  }
}

template <class coeff_t>
inline Vector<coeff_t> operator*(coeff_t alpha, Vector<coeff_t> const &v) {
  Vector<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<std::complex<coeff_t>>
operator*(coeff_t alpha, Vector<std::complex<coeff_t>> const &v) {
  Vector<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(alpha), res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator*(coeff_t alpha, VectorView<coeff_t> v) {
  Vector<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<std::complex<coeff_t>>
operator*(coeff_t alpha, VectorView<std::complex<coeff_t>> v) {
  Vector<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(alpha), res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator*(Vector<coeff_t> const &v, coeff_t alpha) {
  Vector<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<std::complex<coeff_t>>
operator*(Vector<std::complex<coeff_t>> const &v, coeff_t alpha) {
  Vector<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(alpha), res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator*(VectorView<coeff_t> v, coeff_t alpha) {
  Vector<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<std::complex<coeff_t>>
operator*(VectorView<std::complex<coeff_t>> v, coeff_t alpha) {
  Vector<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(alpha), res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator/(Vector<coeff_t> const &v, coeff_t alpha) {
  Vector<coeff_t> res(v);
  Scale(1. / alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<std::complex<coeff_t>>
operator/(Vector<std::complex<coeff_t>> const &v, coeff_t alpha) {
  Vector<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(1. / alpha), res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator/(VectorView<coeff_t> v, coeff_t alpha) {
  Vector<coeff_t> res(v);
  Scale(1. / alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<std::complex<coeff_t>>
operator/(VectorView<std::complex<coeff_t>> v, coeff_t alpha) {
  Vector<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(1. / alpha), res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> &operator*=(Vector<coeff_t> &X, coeff_t alpha) {
  Scale(alpha, X);
  return X;
}

template <class coeff_t>
inline Vector<std::complex<coeff_t>> &
operator*=(Vector<std::complex<coeff_t>> &X, coeff_t alpha) {
  Scale(std::complex<coeff_t>(alpha), X);
  return X;
}

template <class coeff_t>
inline VectorView<coeff_t> operator*=(VectorView<coeff_t> X, coeff_t alpha) {
  Scale(alpha, X);
  return X;
}

template <class coeff_t>
inline VectorView<std::complex<coeff_t>>
operator*=(VectorView<std::complex<coeff_t>> X, coeff_t alpha) {
  Scale(std::complex<coeff_t>(alpha), X);
  return X;
}

template <class coeff_t>
inline Vector<coeff_t> &operator/=(Vector<coeff_t> &X, coeff_t alpha) {
  Scale(1. / alpha, X);
  return X;
}

template <class coeff_t>
inline Vector<std::complex<coeff_t>> &
operator/=(Vector<std::complex<coeff_t>> &X, coeff_t alpha) {
  Scale(std::complex<coeff_t>(1. / alpha, 0), X);
  return X;
}

template <class coeff_t>
inline VectorView<coeff_t> operator/=(VectorView<coeff_t> X, coeff_t alpha) {
  Scale(1. / alpha, X);
  return X;
}

template <class coeff_t>
inline VectorView<std::complex<coeff_t>>
operator/=(VectorView<std::complex<coeff_t>> X, coeff_t alpha) {
  Scale(std::complex<coeff_t>(1. / alpha), X);
  return X;
}

///////

template <class coeff_t>
inline Matrix<coeff_t> operator*(coeff_t alpha, Matrix<coeff_t> const &v) {
  Matrix<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<std::complex<coeff_t>>
operator*(coeff_t alpha, Matrix<std::complex<coeff_t>> const &v) {
  Matrix<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(alpha), res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*(coeff_t alpha, MatrixView<coeff_t> v) {
  Matrix<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<std::complex<coeff_t>>
operator*(coeff_t alpha, MatrixView<std::complex<coeff_t>> v) {
  Matrix<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(alpha), res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*(Matrix<coeff_t> const &v, coeff_t alpha) {
  Matrix<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<std::complex<coeff_t>>
operator*(Matrix<std::complex<coeff_t>> const &v, coeff_t alpha) {
  Matrix<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(alpha), res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*(MatrixView<coeff_t> v, coeff_t alpha) {
  Matrix<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<std::complex<coeff_t>>
operator*(MatrixView<std::complex<coeff_t>> v, coeff_t alpha) {
  Matrix<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(alpha), res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator/(Matrix<coeff_t> const &v, coeff_t alpha) {
  Matrix<coeff_t> res(v);
  Scale(1. / alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<std::complex<coeff_t>>
operator/(Matrix<std::complex<coeff_t>> const &v, coeff_t alpha) {
  Matrix<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(1. / alpha), res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator/(MatrixView<coeff_t> v, coeff_t alpha) {
  Matrix<coeff_t> res(v);
  Scale(1. / alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<std::complex<coeff_t>>
operator/(MatrixView<std::complex<coeff_t>> v, coeff_t alpha) {
  Matrix<std::complex<coeff_t>> res(v);
  Scale(std::complex<coeff_t>(1. / alpha), res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> &operator*=(Matrix<coeff_t> &X, coeff_t alpha) {
  Scale(alpha, X);
  return X;
}

template <class coeff_t>
inline Matrix<std::complex<coeff_t>> &
operator*=(Matrix<std::complex<coeff_t>> &X, coeff_t alpha) {
  Scale(std::complex<coeff_t>(alpha), X);
  return X;
}

template <class coeff_t>
inline MatrixView<coeff_t> operator*=(MatrixView<coeff_t> X, coeff_t alpha) {
  Scale(alpha, X);
  return X;
}

template <class coeff_t>
inline MatrixView<std::complex<coeff_t>>
operator*=(MatrixView<std::complex<coeff_t>> X, coeff_t alpha) {
  Scale(std::complex<coeff_t>(alpha), X);
  return X;
}

template <class coeff_t>
inline Matrix<coeff_t> &operator/=(Matrix<coeff_t> &X, coeff_t alpha) {
  Scale(1. / alpha, X);
  return X;
}

template <class coeff_t>
inline Matrix<std::complex<coeff_t>> &
operator/=(Matrix<std::complex<coeff_t>> &X, coeff_t alpha) {
  Scale(std::complex<coeff_t>(1. / alpha, 0), X);
  return X;
}

template <class coeff_t>
inline MatrixView<coeff_t> operator/=(MatrixView<coeff_t> X, coeff_t alpha) {
  Scale(1. / alpha, X);
  return X;
}

template <class coeff_t>
inline MatrixView<std::complex<coeff_t>>
operator/=(MatrixView<std::complex<coeff_t>> X, coeff_t alpha) {
  Scale(std::complex<coeff_t>(1. / alpha), X);
  return X;
}

} // namespace lila
