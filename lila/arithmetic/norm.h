#pragma once

#include <cassert>

#include <lila/arithmetic/scale.h>
#include <lila/blaslapack/blaslapack.h>

namespace lila {

template <class coeff_t> inline real_t<coeff_t> Norm(Vector<coeff_t> const &v) {
  blas_size_t n = v.n();
  blas_size_t inc = 1;
  return blaslapack::nrm2(&n, LILA_BLAS_CONST_CAST(coeff_t, v.data()), &inc);
}

template <class coeff_t>
inline real_t<coeff_t> Norm(VectorView<coeff_t> const &v) {
  blas_size_t n = v.n();
  blas_size_t inc = v.inc();
  return blaslapack::nrm2(&n, LILA_BLAS_CONST_CAST(coeff_t, v.data()), &inc);
}

template <class coeff_t> inline real_t<coeff_t> Norm(Matrix<coeff_t> const &A) {
  blas_size_t size = A.size();
  blas_size_t inc = 1;
  return blaslapack::nrm2(&size, LILA_BLAS_CONST_CAST(coeff_t, A.data()), &inc);
}

template <class coeff_t>
inline real_t<coeff_t> Norm(MatrixView<coeff_t> const &A) {
  blas_size_t ld = A.ld();
  blas_size_t m = A.m();
  blas_size_t n = A.n();
  blas_size_t incm = A.incm();
  blas_size_t incn = A.incn();

  // Compute norm column-wise
  real_t<coeff_t> norm = 0;
  for (int col = 0; col < n; ++col) {
    real_t<coeff_t> norm_col = blaslapack::nrm2(
        &m, LILA_BLAS_CONST_CAST(coeff_t, A.data() + col * incn * ld), &incm);
    norm += norm_col * norm_col;
  }
  return std::sqrt(norm);
}

template <class coeff_t> inline void Normalize(Vector<coeff_t> &v) {
  coeff_t norm = Norm(v);
  v /= norm;
}

template <class coeff_t> inline void Normalize(VectorView<coeff_t> &v) {
  coeff_t norm = Norm(v);
  v /= norm;
}

template <class coeff_t> inline void Normalize(Matrix<coeff_t> &A) {
  coeff_t norm = Norm(A);
  A /= norm;
}

template <class coeff_t> inline void Normalize(MatrixView<coeff_t> &A) {
  coeff_t norm = Norm(A);
  A /= norm;
}

} // namespace lila
