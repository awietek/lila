#pragma once

#include <lila/blaslapack/blaslapack.h>
#include <lila/matrix.h>

namespace lila {

template <class coeff_t>
inline void CholeskyInplace(Matrix<coeff_t> &A, char uplo = 'U') {
  // check / get dimensions
  assert(A.nrows() == A.ncols());
  blas_size_t n = A.nrows();
  blas_size_t lda = n;
  blas_size_t info = 0;
  blaslapack::potrf(&uplo, &n, LILA_BLAS_CAST(coeff_t, A.data()), &lda, &info);
}

template <class coeff_t>
inline Matrix<coeff_t> Cholesky(const Matrix<coeff_t> &A, char uplo = 'U') {
  auto res = A;
  CholeskyInplace(res, uplo);

  // Set lower triangular to zero
  lila_size_t n = A.nrows();
  if (uplo == 'U') {
    for (lila_size_t j = 0; j < n; ++j)
      for (lila_size_t i = j + 1; i < n; ++i)
        res(i, j) = 0;
  } else if (uplo == 'L') {
    for (lila_size_t i = 0; i < n; ++i)
      for (lila_size_t j = i + 1; j < n; ++j)
        res(i, j) = 0;
  }
  return res;
}

} // namespace lila
