#pragma once

#include <lila/blaslapack/blaslapack.h>
#include <lila/matrix.h>
#include <lila/numeric/complex.h>
#include <lila/special/special.h>
#include <lila/vector.h>

namespace lila {

template <class coeff_t>
inline void Mult(Matrix<coeff_t> const &A, Matrix<coeff_t> const &B,
                 Matrix<coeff_t> &C, coeff_t alpha = 1., coeff_t beta = 0.,
                 char transa = 'N', char transb = 'N') {
  using size_type = blaslapack::blas_size_t;

  // matrix dimensions
  size_type m = (transa == 'N') ? A.nrows() : A.ncols();
  size_type n = (transb == 'N') ? B.ncols() : B.nrows();
  size_type ka = (transa == 'N') ? A.ncols() : A.nrows();
  size_type kb = (transb == 'N') ? B.nrows() : B.ncols();

  assert(ka == kb); // Check if valid multiplication dimensions
  if (((C.m() == 0) && (C.n() == 0)))
    C = Zeros<coeff_t>(m, n);

  if ((C.m() != m) || (C.n() != n))
    C.resize(m, n);

  // leading dimensions
  size_type lda = (transa == 'N') ? m : ka;
  size_type ldb = (transb == 'N') ? ka : n;
  size_type ldc = m;

  blaslapack::gemm(
      &transa, &transb, &m, &n, &ka, LILA_BLAS_CAST(coeff_t, &alpha),
      LILA_BLAS_CONST_CAST(coeff_t, A.data()), &lda,
      LILA_BLAS_CONST_CAST(coeff_t, B.data()), &ldb,
      LILA_BLAS_CAST(coeff_t, &beta), LILA_BLAS_CAST(coeff_t, C.data()), &ldc);
}

template <class coeff_t>
inline Matrix<coeff_t> Mult(Matrix<coeff_t> const &A,
                            Matrix<coeff_t> const &B) {
  Matrix<coeff_t> C;
  Mult(A, B, C);
  return C;
}

template <class coeff_t>
inline void Mult(Matrix<coeff_t> const &A, Vector<coeff_t> const &X,
                 Vector<coeff_t> &Y, coeff_t alpha = 1., coeff_t beta = 0.,
                 char trans = 'N') {
  using size_type = blaslapack::blas_size_t;

  // matrix dimensions
  size_type m = (trans == 'N') ? A.nrows() : A.ncols();
  size_type n = (trans == 'N') ? A.ncols() : A.nrows();

  assert(n == X.size()); // Check if valid multiplication dimensions

  if (Y.size() == 0)
    Y = Zeros<coeff_t>(n);

  // leading dimensions
  size_type lda = A.nrows();
  size_type inc = 1;
  blaslapack::gemv(&trans, &m, &n, LILA_BLAS_CAST(coeff_t, &alpha),
                   LILA_BLAS_CONST_CAST(coeff_t, A.data()), &lda,
                   LILA_BLAS_CONST_CAST(coeff_t, X.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, &beta),
                   LILA_BLAS_CAST(coeff_t, Y.data()), &inc);
}

template <class coeff_t>
inline Vector<coeff_t> Mult(Matrix<coeff_t> const &A,
                            Vector<coeff_t> const &X) {
  Vector<coeff_t> Y;
  Mult(A, X, Y);
  return Y;
}

template <class coeff_t>
inline void Kron(const Matrix<coeff_t> &A, const Matrix<coeff_t> &B,
                 Matrix<coeff_t> &C) {
  using size_type = blaslapack::blas_size_t;
  size_type m = A.nrows();
  size_type n = A.ncols();
  size_type p = B.nrows();
  size_type q = B.ncols();

  size_type mc = m * p;
  size_type nc = n * q;
  C.resize(mc, nc);
  for (int s = 0; s < n; ++s)
    for (int r = 0; r < m; ++r)
      for (int w = 0; w < q; ++w)
        for (int v = 0; v < p; ++v)
          C(p * r + v, q * s + w) = A(r, s) * B(v, w);
}

} // namespace lila
