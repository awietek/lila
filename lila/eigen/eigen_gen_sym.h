#pragma once

#include <utility>

#include <lila/matrix.h>
#include <lila/vector.h>
#include <lila/numeric/complex.h>

namespace lila {

template <class coeff_t>
inline Vector<real_t<coeff_t>>
EigenGenSymDefInplace(Matrix<coeff_t> &A, Matrix<coeff_t> &B,
                      bool do_eigenvectors = true, char uplo = 'U',
                      int itype = 1) {
  using size_type = blaslapack::blas_size_t;

  // check / get dimensions
  assert(A.nrows() == A.ncols());
  assert(B.nrows() == B.ncols());
  assert(A.nrows() == B.nrows());

  char jobz = do_eigenvectors ? 'V' : 'N';
  size_type n = A.nrows();
  size_type lda = n;
  size_type ldb = n;
  int info = 0;

  // get optimal work size
  Vector<real_t<coeff_t>> w(n);
  size_type lwork = -1;
  std::vector<coeff_t> work(1);
  blaslapack::syev(&jobz, &uplo, &n, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                   LILA_BLAS_CAST(real_t<coeff_t>, w.data()),
                   LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);
  assert(info == 0);
  lwork = static_cast<size_type>(real(work[0]));
  work.resize(lwork);

  // Run eigenvalue computation
  blaslapack::sygv(&itype, &jobz, &uplo, &n, LILA_BLAS_CAST(coeff_t, A.data()),
                   &lda, LILA_BLAS_CAST(coeff_t, B.data()), &ldb,
                   LILA_BLAS_CAST(real_t<coeff_t>, w.data()),
                   LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);

  assert(info == 0);
  return w;
}

template <class coeff_t>
inline std::pair<Vector<real_t<coeff_t>>, Matrix<coeff_t>>
EigenGenSymDef(Matrix<coeff_t> const &A, Matrix<coeff_t> const &B,
               char uplo = 'U', int itype = 1) {
  auto eigenvectors = A;
  auto B_copy = B;
  auto eigenvalues =
      EigenGenSymDefInplace(eigenvectors, B_copy, true, uplo, itype);
  return {eigenvalues, eigenvectors};
}

template <class coeff_t>
inline Vector<real_t<coeff_t>>
EigenvaluesGenSymDef(Matrix<coeff_t> const &A, Matrix<coeff_t> const &B,
                     char uplo = 'U', int itype = 1) {
  auto A_copy = A;
  auto B_copy = B;
  return EigenGenSymDefInplace(A_copy, B_copy, false, uplo, itype);
}

} // namespace lila
