#pragma once

#include <tuple>

#include <lila/matrix.h>
#include <lila/vector.h>
#include <lila/numeric/complex.h>

namespace lila {

template <class coeff_t>
inline std::tuple<Vector<complex_t<coeff_t>>, Matrix<complex_t<coeff_t>>,
                  Matrix<complex_t<coeff_t>>>
EigenInplace(Matrix<coeff_t> &A, bool do_right_eigenvectors = true,
             bool do_left_eigenvectors = true) {
  using size_type = blaslapack::blas_size_t;
  using complex_type = complex_t<coeff_t>;
  assert(A.nrows() == A.ncols());

  char jobvl = do_left_eigenvectors ? 'V' : 'N';
  char jobvr = do_right_eigenvectors ? 'V' : 'N';
  size_type n = A.nrows();
  size_type lda = n;

  size_type ldvl = do_left_eigenvectors ? n : 1;
  size_type ldvr = do_right_eigenvectors ? n : 1;

  auto eigenvalues = Vector<complex_type>(n);
  auto left_eigenvectors = Matrix<complex_type>();
  auto right_eigenvectors = Matrix<complex_type>();

  if (do_left_eigenvectors)
    left_eigenvectors.resize(ldvl, n);
  if (do_right_eigenvectors)
    right_eigenvectors.resize(ldvr, n);

  // get optimal work size
  size_type lwork = -1;
  std::vector<coeff_t> work(1);
  int info = 0;
  blaslapack::geev(&jobvl, &jobvr, &n, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                   LILA_BLAS_CAST(complex_type, eigenvalues.data()),
                   LILA_BLAS_CAST(coeff_t, left_eigenvectors.data()), &ldvl,
                   LILA_BLAS_CAST(coeff_t, right_eigenvectors.data()), &ldvr,
                   LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);
  assert(info == 0);
  lwork = static_cast<size_type>(real(work[0]));
  work.resize(lwork);

  // Run eigenvalue computation
  blaslapack::geev(&jobvl, &jobvr, &n, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                   LILA_BLAS_CAST(complex_type, eigenvalues.data()),
                   LILA_BLAS_CAST(coeff_t, left_eigenvectors.data()), &ldvl,
                   LILA_BLAS_CAST(coeff_t, right_eigenvectors.data()), &ldvr,
                   LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);
  assert(info == 0);
  return {eigenvalues, right_eigenvectors, left_eigenvectors};
}

template <class coeff_t>
inline std::tuple<Vector<complex_t<coeff_t>>, Matrix<complex_t<coeff_t>>,
                  Matrix<complex_t<coeff_t>>>
Eigen(Matrix<coeff_t> const &A, bool do_right_eigenvectors = true,
      bool do_left_eigenvectors = true) {
  auto A_copy = A;
  return EigenInplace(A_copy, do_right_eigenvectors, do_left_eigenvectors);
}

template <class coeff_t>
inline Vector<complex_t<coeff_t>> Eigenvalues(Matrix<coeff_t> const &A) {
  [[maybe_unused]] auto [eigs, l, r] = Eigen(A, false, false);
  return eigs;
}

} // namespace lila
