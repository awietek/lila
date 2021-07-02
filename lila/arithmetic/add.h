#pragma once

#include <cassert>

#include <lila/matrix.h>
#include <lila/vector.h>

#include <lila/blaslapack/blaslapack.h>

namespace lila {

template <class coeff_t>
inline void Add(Vector<coeff_t> const &v, Vector<coeff_t> &w,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  using size_type = blaslapack::blas_size_t;
  assert(v.n() == w.n());
  size_type n = v.n();
  size_type inc = 1;
  blaslapack::axpy(&n, LILA_BLAS_CAST(coeff_t, &alpha),
                   LILA_BLAS_CONST_CAST(coeff_t, v.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, w.data()), &inc);
}

template <class coeff_t>
inline void Add(VectorView<coeff_t> v, VectorView<coeff_t> w,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  using size_type = blaslapack::blas_size_t;
  assert(v.n() == w.n());
  size_type n = v.n();
  size_type v_inc = v.inc();
  size_type w_inc = w.inc();
  blaslapack::axpy(&n, LILA_BLAS_CAST(coeff_t, &alpha),
                   LILA_BLAS_CONST_CAST(coeff_t, v.data()), &v_inc,
                   LILA_BLAS_CAST(coeff_t, w.data()), &w_inc);
}

template <class coeff_t>
inline void Add(Matrix<coeff_t> const &A, Matrix<coeff_t> &B,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  using size_type = blaslapack::blas_size_t;
  assert(A.m() == B.m());
  assert(A.n() == B.n());
  size_type size = A.size();
  size_type inc = 1;
  blaslapack::axpy(&size, LILA_BLAS_CAST(coeff_t, &alpha),
                   LILA_BLAS_CONST_CAST(coeff_t, A.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, B.data()), &inc);
}

template <class coeff_t>
inline void Add(MatrixView<coeff_t> A, MatrixView<coeff_t> B,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  using size_type = blaslapack::blas_size_t;
  assert(A.m() == B.m());
  assert(A.n() == B.n());

  size_type A_incm = A.incm();
  size_type A_incn = A.incn();
  size_type B_incm = B.incm();
  size_type B_incn = B.incn();

  size_type m = A.m();
  size_type n = A.n();
  size_type A_ld = A.ld();
  size_type B_ld = B.ld();

  // Perform a column-wise axpy
  for (int col = 0; col < n; ++col) {
    blaslapack::axpy(
        &m, LILA_BLAS_CAST(coeff_t, &alpha),
        LILA_BLAS_CONST_CAST(coeff_t, A.data() + col * A_incn * A_ld), &A_incm,
        LILA_BLAS_CAST(coeff_t, B.data() + col * B_incn * B_ld), &B_incm);
  }
}

template <class coeff_t>
inline Vector<coeff_t> operator+(Vector<coeff_t> const &v,
                                 Vector<coeff_t> const &w) {
  Vector<coeff_t> res(w); // Creates a new Vector
  Add(v, res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator+(VectorView<coeff_t> v, VectorView<coeff_t> w) {
  Vector<coeff_t> res(w); // Creates a new Vector
  Add(v, res);
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator-(Vector<coeff_t> const &v,
                                 Vector<coeff_t> const &w) {
  Vector<coeff_t> res(v); // Creates a new Vector
  Add(w, res, static_cast<coeff_t>(-1.));
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator-(VectorView<coeff_t> v, VectorView<coeff_t> w) {
  Vector<coeff_t> res(v); // Creates a new Vector
  Add(w, res, static_cast<coeff_t>(-1.));
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator+(Vector<coeff_t> const &v, coeff_t c) {
  Vector<coeff_t> res(v); // Creates a new Vector
  Map(res, [&c](coeff_t &x) { x = x + c; });
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator+(VectorView<coeff_t> v, coeff_t c) {
  Vector<coeff_t> res(v); // Creates a new Vector
  Map(res, [&c](coeff_t &x) { x = x + c; });
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator-(Vector<coeff_t> const &v, coeff_t c) {
  Vector<coeff_t> res(v); // Creates a new Vector
  Map(res, [&c](coeff_t &x) { x = x - c; });
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator-(VectorView<coeff_t> v, coeff_t c) {
  Vector<coeff_t> res(v); // Creates a new Vector
  Map(res, [&c](coeff_t &x) { x = x - c; });
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> &operator+=(Vector<coeff_t> &v,
                                   Vector<coeff_t> const &w) {
  Add(w, v);
  return v;
}

template <class coeff_t>
inline VectorView<coeff_t> operator+=(VectorView<coeff_t> v,
                                      VectorView<coeff_t> w) {
  Add(w, v);
  return v;
}

template <class coeff_t>
inline Vector<coeff_t> &operator+=(Vector<coeff_t> &v, coeff_t c) {
  Map(v, [&c](coeff_t &x) { x = x + c; });
  return v;
}

template <class coeff_t>
inline VectorView<coeff_t> operator+=(VectorView<coeff_t> v, coeff_t c) {
  Map(v, [&c](coeff_t &x) { x = x + c; });
  return v;
}

template <class coeff_t>
inline Vector<coeff_t> &operator-=(Vector<coeff_t> &v,
                                   Vector<coeff_t> const &w) {
  Add(w, v, static_cast<coeff_t>(-1.));
  return v;
}

template <class coeff_t>
inline VectorView<coeff_t> operator-=(VectorView<coeff_t> v,
                                      VectorView<coeff_t> w) {
  Add(w, v, static_cast<coeff_t>(-1.));
  return v;
}

template <class coeff_t>
inline Vector<coeff_t> &operator-=(Vector<coeff_t> &v, coeff_t c) {
  Map(v, [&c](coeff_t &x) { x = x - c; });
  return v;
}

template <class coeff_t>
inline VectorView<coeff_t> operator-=(VectorView<coeff_t> v, coeff_t c) {
  Map(v, [&c](coeff_t &x) { x = x - c; });
  return v;
}

template <class coeff_t>
inline Vector<coeff_t> operator-(Vector<coeff_t> const &v) {
  Vector<coeff_t> res(v);
  Map(res, [](coeff_t &x) { x = -x; });
  return res;
}

template <class coeff_t>
inline Vector<coeff_t> operator-(VectorView<coeff_t> v) {
  Vector<coeff_t> res(v);
  Map(res, [](coeff_t &x) { x = -x; });
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator+(Matrix<coeff_t> const &v,
                                 Matrix<coeff_t> const &w) {
  Matrix<coeff_t> res(w); // Creates a new Matrix
  Add(v, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator+(MatrixView<coeff_t> v, MatrixView<coeff_t> w) {
  Matrix<coeff_t> res(w); // Creates a new Matrix
  Add(v, res);
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator-(Matrix<coeff_t> const &v,
                                 Matrix<coeff_t> const &w) {
  Matrix<coeff_t> res(v); // Creates a new Matrix
  Add(w, res, static_cast<coeff_t>(-1.));
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator-(MatrixView<coeff_t> v, MatrixView<coeff_t> w) {
  Matrix<coeff_t> res(v); // Creates a new Matrix
  Add(w, res, static_cast<coeff_t>(-1.));
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator+(Matrix<coeff_t> const &A, coeff_t c) {
  Matrix<coeff_t> res(A); // Creates a new Matrix
  Map(res, [&c](coeff_t &x) { x = x + c; });
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator+(MatrixView<coeff_t> A, coeff_t c) {
  Matrix<coeff_t> res(A); // Creates a new Matrix
  Map(res, [&c](coeff_t &x) { x = x + c; });
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator-(Matrix<coeff_t> const &A, coeff_t c) {
  Matrix<coeff_t> res(A); // Creates a new Matrix
  Map(res, [&c](coeff_t &x) { x = x - c; });
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator-(MatrixView<coeff_t> A, coeff_t c) {
  Matrix<coeff_t> res(A); // Creates a new Matrix
  Map(res, [&c](coeff_t &x) { x = x - c; });
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> &operator+=(Matrix<coeff_t> &A,
                                   Matrix<coeff_t> const &B) {
  Add(B, A);
  return A;
}

template <class coeff_t>
inline MatrixView<coeff_t> operator+=(MatrixView<coeff_t> A,
                                      MatrixView<coeff_t> B) {
  Add(B, A);
  return A;
}
template <class coeff_t>
inline Matrix<coeff_t> &operator+=(Matrix<coeff_t> &v, coeff_t c) {
  Map(v, [&c](coeff_t &x) { x = x + c; });
  return v;
}

template <class coeff_t>
inline MatrixView<coeff_t> operator+=(MatrixView<coeff_t> v, coeff_t c) {
  Map(v, [&c](coeff_t &x) { x = x + c; });
  return v;
}

template <class coeff_t>
inline Matrix<coeff_t> &operator-=(Matrix<coeff_t> &A,
                                   Matrix<coeff_t> const &B) {
  Add(B, A, static_cast<coeff_t>(-1.));
  return A;
}

template <class coeff_t>
inline MatrixView<coeff_t> operator-=(MatrixView<coeff_t> A,
                                      MatrixView<coeff_t> B) {
  Add(B, A, static_cast<coeff_t>(-1.));
  return A;
}

template <class coeff_t>
inline Matrix<coeff_t> &operator-=(Matrix<coeff_t> &v, coeff_t c) {
  Map(v, [&c](coeff_t &x) { x = x - c; });
  return v;
}

template <class coeff_t>
inline MatrixView<coeff_t> operator-=(MatrixView<coeff_t> v, coeff_t c) {
  Map(v, [&c](coeff_t &x) { x = x - c; });
  return v;
}

template <class coeff_t>
inline Matrix<coeff_t> operator-(Matrix<coeff_t> const &v) {
  Matrix<coeff_t> res(v);
  Map(res, [](coeff_t &x) { x = -x; });
  return res;
}

template <class coeff_t>
inline Matrix<coeff_t> operator-(MatrixView<coeff_t> v) {
  Matrix<coeff_t> res(v);
  Map(res, [](coeff_t &x) { x = -x; });
  return res;
}

} // namespace lila
