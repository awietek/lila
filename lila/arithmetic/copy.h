#pragma once

#include <cassert>

#include <lila/blaslapack/blaslapack.h>
#include <lila/matrix.h>
#include <lila/vector.h>
#include <lila/views/matrix_view.h>
#include <lila/views/vector_view.h>

namespace lila {

template <class coeff_t>
inline void Copy(Vector<coeff_t> const &v, Vector<coeff_t> &w) {
  using size_type = blaslapack::blas_size_t;
  assert(v.n() == w.n());
  size_type n = v.n();
  size_type inc = 1;
  blaslapack::copy(&n, LILA_BLAS_CONST_CAST(coeff_t, v.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, w.data()), &inc);
}

template <class coeff_t>
inline void Copy(VectorView<coeff_t> v, VectorView<coeff_t> w) {
  using size_type = blaslapack::blas_size_t;
  assert(v.n() == w.n());
  size_type n = v.n();
  size_type v_inc = v.inc();
  size_type w_inc = w.inc();
  blaslapack::copy(&n, LILA_BLAS_CONST_CAST(coeff_t, v.data()), &v_inc,
                   LILA_BLAS_CAST(coeff_t, w.data()), &w_inc);
}

template <class coeff_t>
inline void Copy(Matrix<coeff_t> const &A, Matrix<coeff_t> &B) {
  using size_type = blaslapack::blas_size_t;
  assert(A.m() == B.m());
  assert(A.n() == B.n());
  size_type size = A.size();
  size_type inc = 1;
  blaslapack::copy(&size, LILA_BLAS_CONST_CAST(coeff_t, A.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, B.data()), &inc);
}

template <class coeff_t>
inline void Copy(MatrixView<coeff_t> A, MatrixView<coeff_t> B) {
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

  // Perform a column-wise copy
  for (int col = 0; col < n; ++col) {
    blaslapack::copy(
        &m, LILA_BLAS_CONST_CAST(coeff_t, A.data() + col * A_incn * A_ld),
        &A_incm, LILA_BLAS_CAST(coeff_t, B.data() + col * B_incn * B_ld),
        &B_incm);
  }
}

} // namespace lila
