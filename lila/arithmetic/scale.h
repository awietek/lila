#pragma once

#include <cassert>

#include <lila/matrix.h>
#include <lila/vector.h>

#include <lila/blaslapack/blaslapack.h>

namespace lila {

template <class coeff_t> inline void Scale(coeff_t alpha, Vector<coeff_t> &v) {
  using size_type = blaslapack::blas_size_t;
  size_type size = v.size();
  size_type inc = 1;
  blaslapack::scal(&size, LILA_BLAS_CONST_CAST(coeff_t, &alpha),
                   LILA_BLAS_CAST(coeff_t, v.data()), &inc);
}

template <class coeff_t>
inline void Scale(coeff_t alpha, VectorView<coeff_t> v) {
  using size_type = blaslapack::blas_size_t;
  size_type size = v.size();
  size_type inc = v.inc();
  blaslapack::scal(&size, LILA_BLAS_CONST_CAST(coeff_t, &alpha),
                   LILA_BLAS_CAST(coeff_t, v.data()), &inc);
}

template <class coeff_t> inline void Scale(coeff_t alpha, Matrix<coeff_t> &A) {
  using size_type = blaslapack::blas_size_t;
  size_type dx = A.size();
  size_type inc = 1;
  blaslapack::scal(&dx, LILA_BLAS_CONST_CAST(coeff_t, &alpha),
                   LILA_BLAS_CAST(coeff_t, A.data()), &inc);
}

template <class coeff_t>
inline void Scale(coeff_t alpha, MatrixView<coeff_t> A) {
  using size_type = blaslapack::blas_size_t;

  size_type ld = A.ld();
  size_type m = A.m();
  size_type n = A.n();
  size_type incm = A.incm();
  size_type incn = A.incn();

  // Perform a column-wise scale
  for (int col = 0; col < n; ++col) {
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
inline Vector<coeff_t> operator*(coeff_t alpha, VectorView<coeff_t> v) {
  Vector<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator*(Vector<coeff_t> const &v, coeff_t alpha) {
  Vector<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator*(VectorView<coeff_t> v, coeff_t alpha) {
  Vector<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator/(Vector<coeff_t> const &v, coeff_t alpha) {
  Vector<coeff_t> res(v);
  Scale(1. / alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator/(VectorView<coeff_t> v, coeff_t alpha) {
  Vector<coeff_t> res(v);
  Scale(1. / alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> &operator*=(Vector<coeff_t> &X, coeff_t alpha) {
  Scale(alpha, X);
  return X;
}

template <class coeff_t>
inline VectorView<coeff_t> operator*=(VectorView<coeff_t> X, coeff_t alpha) {
  Scale(alpha, X);
  return X;
}

template <class coeff_t>
inline Vector<coeff_t> &operator/=(Vector<coeff_t> &X, coeff_t alpha) {
  Scale(1. / alpha, X);
  return X;
}

template <class coeff_t>
inline VectorView<coeff_t> operator/=(VectorView<coeff_t> X, coeff_t alpha) {
  Scale(1. / alpha, X);
  return X;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*(coeff_t alpha, Matrix<coeff_t> const &v) {
  Matrix<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*(coeff_t alpha, MatrixView<coeff_t> v) {
  Matrix<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*(Matrix<coeff_t> const &v, coeff_t alpha) {
  Matrix<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*(MatrixView<coeff_t> v, coeff_t alpha) {
  Matrix<coeff_t> res(v);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator/(Matrix<coeff_t> const &v, coeff_t alpha) {
  Matrix<coeff_t> res(v);
  Scale(1. / alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator/(MatrixView<coeff_t> v, coeff_t alpha) {
  Matrix<coeff_t> res(v);
  Scale(1. / alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> &operator*=(Matrix<coeff_t> &X, coeff_t alpha) {
  Scale(alpha, X);
  return X;
}

template <class coeff_t>
inline MatrixView<coeff_t> operator*=(MatrixView<coeff_t> X, coeff_t alpha) {
  Scale(alpha, X);
  return X;
}

template <class coeff_t>
inline Matrix<coeff_t> &operator/=(Matrix<coeff_t> &X, coeff_t alpha) {
  Scale(1. / alpha, X);
  return X;
}

template <class coeff_t>
inline MatrixView<coeff_t> operator/=(MatrixView<coeff_t> X, coeff_t alpha) {
  Scale(1. / alpha, X);
  return X;
}

} // namespace lila
