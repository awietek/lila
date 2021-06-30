#pragma once

#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>

#include <lila/matrix.h>
#include <lila/numeric/compare.h>
#include <lila/numeric/complex.h>
#include <lila/vector.h>

#include <lila/blaslapack/blaslapack.h>

namespace lila {

template <class coeff_t>
inline void Add(const Vector<coeff_t> &X, Vector<coeff_t> &Y,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  using size_type = blaslapack::blas_size_t;
  size_type dx = X.size();
  size_type dy = Y.size();
  assert(dx == dy); // Check if valid dimensions
  size_type inc = 1;
  blaslapack::axpy(&dx, LILA_BLAS_CAST(coeff_t, &alpha),
                   LILA_BLAS_CONST_CAST(coeff_t, X.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, Y.data()), &inc);
}

template <class coeff_t>
inline void Add(const Matrix<coeff_t> &X, Matrix<coeff_t> &Y,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  using size_type = blaslapack::blas_size_t;
  size_type dx = X.size();
  size_type dy = Y.size();
  assert(dx == dy); // Check if valid dimensions
  size_type inc = 1;
  blaslapack::axpy(&dx, LILA_BLAS_CAST(coeff_t, &alpha),
                   LILA_BLAS_CONST_CAST(coeff_t, X.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, Y.data()), &inc);
}

template <class coeff_t>
inline void Scale(const coeff_t &alpha, Vector<coeff_t> &X) {
  using size_type = blaslapack::blas_size_t;
  size_type dx = X.size();
  size_type inc = 1;
  blaslapack::scal(&dx, LILA_BLAS_CONST_CAST(coeff_t, &alpha),
                   LILA_BLAS_CAST(coeff_t, X.data()), &inc);
}

template <class coeff_t>
inline void Scale(const coeff_t &alpha, Matrix<coeff_t> &X) {
  using size_type = blaslapack::blas_size_t;
  size_type dx = X.size();
  size_type inc = 1;
  blaslapack::scal(&dx, LILA_BLAS_CONST_CAST(coeff_t, &alpha),
                   LILA_BLAS_CAST(coeff_t, X.data()), &inc);
}

template <class coeff_t>
inline coeff_t Dot(const Vector<coeff_t> &X, const Vector<coeff_t> &Y) {
  using size_type = blaslapack::blas_size_t;
  size_type dx = X.size();
  size_type dy = Y.size();
  assert(dx == dy); // Check if valid dimensions
  size_type inc = 1;
  blaslapack::blas_t<coeff_t> result =
      blaslapack::dot(&dx, LILA_BLAS_CONST_CAST(coeff_t, X.data()), &inc,
                      LILA_BLAS_CONST_CAST(coeff_t, Y.data()), &inc);
  return blaslapack::blas_to_lila(result);
}

template <class coeff_t>
inline coeff_t Dot(const Matrix<coeff_t> &X, const Matrix<coeff_t> &Y) {
  using size_type = blaslapack::blas_size_t;
  size_type dx = X.size();
  size_type dy = Y.size();
  assert(dx == dy); // Check if valid dimensions
  size_type inc = 1;
  blaslapack::blas_t<coeff_t> result =
      blaslapack::dot(&dx, LILA_BLAS_CONST_CAST(coeff_t, X.data()), &inc,
                      LILA_BLAS_CONST_CAST(coeff_t, Y.data()), &inc);
  return blaslapack::blas_to_lila(result);
}

template <class coeff_t> inline real_t<coeff_t> Norm(const Vector<coeff_t> &X) {
  // TODO: use proper LAPACK function here
  return sqrt(real(Dot(X, X)));
}

template <class coeff_t> inline real_t<coeff_t> Norm(const Matrix<coeff_t> &X) {
  // TODO: use proper LAPACK function here
  return sqrt(real(Dot(X, X)));
}

template <class coeff_t> inline void Normalize(Vector<coeff_t> &X) {
  coeff_t nrm = Norm(X);
  Scale((coeff_t)1. / nrm, X);
}

template <class coeff_t> inline void Normalize(Matrix<coeff_t> &X) {
  coeff_t nrm = Norm(X);
  Scale((coeff_t)1. / nrm, X);
}

template <class coeff_t>
inline real_t<coeff_t> NormLi(const Matrix<coeff_t> &X) {
  real_t<coeff_t> value = 0.0;
  for (int i = 0; i < X.nrows(); i++) {
    real_t<coeff_t> row_sum = 0.0;
    for (int j = 0; j < X.ncols(); j++) {
      row_sum += std::abs(X(i, j));
    }
    value = std::max(value, row_sum);
  }
  return value;
}

template <class coeff_t> inline real_t<coeff_t> log2abs(const coeff_t &x) {
  real_t<coeff_t> value;

  if (x == 0.0) {
    value = -1e30;
  } else {
    value = log(std::abs(x)) / log(2.0);
  }

  return value;
}

template <class coeff_t> inline coeff_t Sum(const Vector<coeff_t> &X) {
  return std::accumulate(X.begin(), X.end(), 0.);
}
template <class coeff_t> inline coeff_t Sum(const Matrix<coeff_t> &X) {
  return std::accumulate(X.begin(), X.end(), 0.);
}

template <class coeff_t>
inline Vector<coeff_t> operator+(const Vector<coeff_t> &X,
                                 const Vector<coeff_t> &Y) {
  Vector<coeff_t> res(Y);
  Add(X, res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator-(const Vector<coeff_t> &X,
                                 const Vector<coeff_t> &Y) {
  Vector<coeff_t> res(X);
  Add(Y, res, static_cast<coeff_t>(-1.));
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator+(const Vector<coeff_t> &X, const coeff_t &c) {
  Vector<coeff_t> res(X);
  Map(res, [&c](coeff_t &x) { x = x + c; });
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> &operator+=(Vector<coeff_t> &a,
                                   const Vector<coeff_t> &b) {
  a = a + b;
  return a;
}
template <class coeff_t>
inline Vector<coeff_t> &operator-=(Vector<coeff_t> &a,
                                   const Vector<coeff_t> &b) {
  a = a - b;
  return a;
}

template <class coeff_t>
inline Vector<coeff_t> operator-(const Vector<coeff_t> &X, const coeff_t &c) {
  Vector<coeff_t> res(X);
  Map(res, [&c](coeff_t &x) { x = x - c; });
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator-(const Vector<coeff_t> &X) {
  Vector<coeff_t> res(X);
  Scale(static_cast<coeff_t>(-1.), res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator*(const coeff_t &alpha,
                                 const Vector<coeff_t> &X) {
  Vector<coeff_t> res(X);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator*=(Vector<coeff_t> &X, const coeff_t &alpha) {
  X = alpha * X;
  return X;
}

template <class coeff_t>
inline Vector<coeff_t> operator/=(Vector<coeff_t> &X, const coeff_t &alpha) {
  X = X / alpha;
  return X;
}

template <class coeff_t>
inline Vector<coeff_t> operator*(const Vector<coeff_t> &X,
                                 const coeff_t &alpha) {
  return operator*(alpha, X);
}

template <class coeff_t>
inline Vector<coeff_t> operator/(const Vector<coeff_t> &X,
                                 const coeff_t &alpha) {
  assert(!close(alpha, static_cast<coeff_t>(0.)));
  coeff_t invalpha = static_cast<coeff_t>(1.) / alpha;
  return operator*(invalpha, X);
}

template <class coeff_t>
inline Matrix<coeff_t> operator+(const Matrix<coeff_t> &X,
                                 const Matrix<coeff_t> &Y) {
  Matrix<coeff_t> res(Y);
  Add(X, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator-(const Matrix<coeff_t> &X,
                                 const Matrix<coeff_t> &Y) {
  Matrix<coeff_t> res(X);
  Add(Y, res, static_cast<coeff_t>(-1.));
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator+(const Matrix<coeff_t> &X, const coeff_t &c) {
  Matrix<coeff_t> res(X);
  Map(res, [&c](coeff_t &x) { x = x + c; });
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> &operator+=(Matrix<coeff_t> &a,
                                   const Matrix<coeff_t> &b) {
  a = a + b;
  return a;
}
template <class coeff_t>
inline Matrix<coeff_t> &operator-=(Matrix<coeff_t> &a,
                                   const Matrix<coeff_t> &b) {
  a = a - b;
  return a;
}

template <class coeff_t>
inline Matrix<coeff_t> operator-(const Matrix<coeff_t> &X, const coeff_t &c) {
  Matrix<coeff_t> res(X);
  Map(res, [&c](coeff_t &x) { x = x - c; });
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator-(const Matrix<coeff_t> &X) {
  Matrix<coeff_t> res(X);
  Scale(static_cast<coeff_t>(-1.), res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*(const coeff_t &alpha,
                                 const Matrix<coeff_t> &X) {
  Matrix<coeff_t> res(X);
  Scale(alpha, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*=(Matrix<coeff_t> &X, const coeff_t &alpha) {
  X = alpha * X;
  return X;
}

template <class coeff_t>
inline Matrix<coeff_t> operator/=(Matrix<coeff_t> &X, const coeff_t &alpha) {
  X = X / alpha;
  return X;
}

template <class coeff_t>
inline Matrix<coeff_t> operator*(const Matrix<coeff_t> &X,
                                 const coeff_t &alpha) {
  return operator*(alpha, X);
}

template <class coeff_t>
inline Matrix<coeff_t> operator/(const Matrix<coeff_t> &X,
                                 const coeff_t &alpha) {
  assert(!close(alpha, static_cast<coeff_t>(0.)));
  coeff_t invalpha = static_cast<coeff_t>(1.) / alpha;
  return operator*(invalpha, X);
}

} // namespace lila
