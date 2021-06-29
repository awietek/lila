#pragma once

#include <utility>

#include <lila/matrix.h>
#include <lila/vector.h>
#include <lila/numeric/complex.h>

namespace lila {

template <class coeff_t>
inline Vector<real_t<coeff_t>> EigenSymInplace(Matrix<coeff_t> &A,
                                               bool do_eigenvectors = true,
                                               char uplo = 'U') {
  using size_type = blaslapack::blas_size_t;
  assert(A.nrows() == A.ncols());

  char jobz = do_eigenvectors ? 'V' : 'N';
  size_type n = A.nrows();
  size_type lda = n;

  // get optimal work size
  Vector<real_t<coeff_t>> w(n);
  size_type lwork = -1;
  std::vector<coeff_t> work(1);
  int info = 0;

  blaslapack::syev(&jobz, &uplo, &n, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                   LILA_BLAS_CAST(real_t<coeff_t>, w.data()),
                   LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);

  assert(info == 0);
  lwork = static_cast<size_type>(real(work[0]));
  work.resize(lwork);

  // Run eigenvalue computation
  blaslapack::syev(&jobz, &uplo, &n, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                   LILA_BLAS_CAST(real_t<coeff_t>, w.data()),
                   LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);

  assert(info == 0);
  return w;
}

template <class coeff_t>
inline std::pair<Vector<real_t<coeff_t>>, Matrix<coeff_t>>
EigenSym(Matrix<coeff_t> const &A, char uplo = 'U') {
  auto eigenvectors = A;
  auto eigenvalues = EigenSymInplace(eigenvectors, true, uplo);
  return {eigenvalues, eigenvectors};
}

template <class coeff_t>
inline Vector<real_t<coeff_t>> EigenvaluesSym(Matrix<coeff_t> const &A,
                                              char uplo = 'U') {
  Matrix<coeff_t> A_copy = A;
  return EigenSymInplace(A_copy, false, uplo);
}

} // namespace lila
