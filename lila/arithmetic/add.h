#pragma once

#include <cassert>

#include <lila/matrix.h>
#include <lila/vector.h>
#include <lila/views/matrix_view.h>
#include <lila/views/vector_view.h>

#include <lila/blaslapack/blaslapack.h>

namespace lila {

template <class coeff_t>
inline void Add(Vector<coeff_t> const &v, Vector<coeff_t> &w,
                coeff_t alpha = static_cast<coeff_t>(1.)) {

  assert(v.n() == w.n());
  blas_size_t n = v.n();
  blas_size_t inc = 1;
  blaslapack::axpy(&n, LILA_BLAS_CAST(coeff_t, &alpha),
                   LILA_BLAS_CONST_CAST(coeff_t, v.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, w.data()), &inc);
}

template <class coeff_t>
inline void Add(VectorView<coeff_t> v, VectorView<coeff_t> w,
                coeff_t alpha = static_cast<coeff_t>(1.)) {

  assert(v.n() == w.n());
  blas_size_t n = v.n();
  blas_size_t v_inc = v.inc();
  blas_size_t w_inc = w.inc();
  blaslapack::axpy(&n, LILA_BLAS_CAST(coeff_t, &alpha),
                   LILA_BLAS_CONST_CAST(coeff_t, v.data()), &v_inc,
                   LILA_BLAS_CAST(coeff_t, w.data()), &w_inc);
}

template <class coeff_t>
inline void Add(Vector<coeff_t> const &v, VectorView<coeff_t> w,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  Add(VectorView<coeff_t>(v), w, alpha);
}

template <class coeff_t>
inline void Add(VectorView<coeff_t> v, Vector<coeff_t> &w,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  Add(v, VectorView<coeff_t>(w), alpha);
}

template <class coeff_t>
inline void Add(Matrix<coeff_t> const &A, Matrix<coeff_t> &B,
                coeff_t alpha = static_cast<coeff_t>(1.)) {

  assert(A.m() == B.m());
  assert(A.n() == B.n());
  blas_size_t size = A.size();
  blas_size_t inc = 1;
  blaslapack::axpy(&size, LILA_BLAS_CAST(coeff_t, &alpha),
                   LILA_BLAS_CONST_CAST(coeff_t, A.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, B.data()), &inc);
}

template <class coeff_t>
inline void Add(MatrixView<coeff_t> A, MatrixView<coeff_t> B,
                coeff_t alpha = static_cast<coeff_t>(1.)) {

  assert(A.m() == B.m());
  assert(A.n() == B.n());

  blas_size_t A_incm = A.incm();
  blas_size_t A_incn = A.incn();
  blas_size_t B_incm = B.incm();
  blas_size_t B_incn = B.incn();

  blas_size_t m = A.m();
  blas_size_t n = A.n();
  blas_size_t A_ld = A.ld();
  blas_size_t B_ld = B.ld();

  // Perform a column-wise axpy
  for (blas_size_t col = 0; col < n; ++col) {
    blaslapack::axpy(
        &m, LILA_BLAS_CAST(coeff_t, &alpha),
        LILA_BLAS_CONST_CAST(coeff_t, A.data() + col * A_incn * A_ld), &A_incm,
        LILA_BLAS_CAST(coeff_t, B.data() + col * B_incn * B_ld), &B_incm);
  }
}

template <class coeff_t>
inline void Add(Matrix<coeff_t> const &A, MatrixView<coeff_t> B,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  Add(MatrixView<coeff_t>(A), B, alpha);
}

template <class coeff_t>
inline void Add(MatrixView<coeff_t> A, Matrix<coeff_t> &B,
                coeff_t alpha = static_cast<coeff_t>(1.)) {
  Add(A, MatrixView<coeff_t>(B), alpha);
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
