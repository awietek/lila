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

  assert(v.n() == w.n());
  blas_size_t n = v.n();
  blas_size_t inc = 1;
  blaslapack::copy(&n, LILA_BLAS_CONST_CAST(coeff_t, v.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, w.data()), &inc);
}

template <class coeff_t>
inline void Copy(VectorView<coeff_t> v, VectorView<coeff_t> w) {

  assert(v.n() == w.n());
  blas_size_t n = v.n();
  blas_size_t v_inc = v.inc();
  blas_size_t w_inc = w.inc();
  blaslapack::copy(&n, LILA_BLAS_CONST_CAST(coeff_t, v.data()), &v_inc,
                   LILA_BLAS_CAST(coeff_t, w.data()), &w_inc);
}

template <class coeff_t>
inline void Copy(VectorView<coeff_t> v, Vector<coeff_t> &w) {
  Copy(v, VectorView<coeff_t>(w));
}

template <class coeff_t>
inline void Copy(Vector<coeff_t> const &v, VectorView<coeff_t> w) {
  Copy(VectorView<coeff_t>(v), w);
}

template <class coeff_t>
inline void Copy(Matrix<coeff_t> const &A, Matrix<coeff_t> &B) {

  assert(A.m() == B.m());
  assert(A.n() == B.n());
  blas_size_t size = A.size();
  blas_size_t inc = 1;
  blaslapack::copy(&size, LILA_BLAS_CONST_CAST(coeff_t, A.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, B.data()), &inc);
}

template <class coeff_t>
inline void Copy(MatrixView<coeff_t> A, MatrixView<coeff_t> B) {

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

  // Perform a column-wise copy
  for (blas_size_t col = 0; col < n; ++col) {
    blaslapack::copy(
        &m, LILA_BLAS_CONST_CAST(coeff_t, A.data() + col * A_incn * A_ld),
        &A_incm, LILA_BLAS_CAST(coeff_t, B.data() + col * B_incn * B_ld),
        &B_incm);
  }
}

template <class coeff_t>
inline void Copy(MatrixView<coeff_t> v, Matrix<coeff_t> &w) {
  Copy(v, MatrixView<coeff_t>(w));
}

template <class coeff_t>
inline void Copy(Matrix<coeff_t> const &v, MatrixView<coeff_t> w) {
  Copy(MatrixView<coeff_t>(v), w);
}

} // namespace lila
