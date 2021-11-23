#pragma once

#include <lila/blaslapack/blaslapack.h>
#include <lila/matrix.h>
#include <lila/special/special.h>
#include <lila/vector.h>

namespace lila {

template <class coeff_t>
inline std::vector<blas_size_t> Solve(Matrix<coeff_t> &A, Matrix<coeff_t> &B) {
  // check / get dimensions
  assert(A.nrows() == A.ncols());
  assert(A.nrows() == B.nrows());
  blas_size_t n = A.nrows();
  blas_size_t n_rhs = B.ncols();
  blas_size_t lda = n;
  blas_size_t ldb = n;

  std::vector<blas_size_t> ipiv(n);
  blas_size_t info = 0;

  blaslapack::gesv(&n, &n_rhs, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                   ipiv.data(), LILA_BLAS_CAST(coeff_t, B.data()), &ldb, &info);

  return ipiv;
}

template <class coeff_t>
inline std::vector<blas_size_t> Solve(Matrix<coeff_t> &A, Vector<coeff_t> &X) {

  // check / get dimensions
  assert(A.nrows() == X.n());
  blas_size_t n = A.nrows();
  blas_size_t n_rhs = 1;
  blas_size_t lda = n;
  blas_size_t ldb = n;

  std::vector<blas_size_t> ipiv(n);
  blas_size_t info = 0;

  blaslapack::gesv(&n, &n_rhs, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                   ipiv.data(), LILA_BLAS_CAST(coeff_t, X.data()), &ldb, &info);
  assert(info == 0);

  return ipiv;
}

template <class coeff_t>
inline std::vector<blas_size_t> LUDecompose(Matrix<coeff_t> &A) {

  // check / get dimensions
  blas_size_t m = A.nrows();
  blas_size_t n = A.ncols();
  blas_size_t lda = m;

  std::vector<blas_size_t> ipiv(n);
  blas_size_t info = 0;

  blaslapack::getrf(&m, &n, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                    ipiv.data(), &info);
  // assert(info == 0);
  // if (info != 0)
  //   printf("info %d\n", info);

  return ipiv;
}

template <class coeff_t>
inline void LUSolve(const Matrix<coeff_t> &A,
                    const std::vector<blas_size_t> &ipiv, Matrix<coeff_t> &B,
                    char trans = 'N') {
  // check / get dimensions
  assert(A.nrows() == A.ncols());
  assert(A.nrows() == B.nrows());
  blas_size_t n = A.nrows();
  blas_size_t n_rhs = B.ncols();
  blas_size_t lda = n;
  blas_size_t ldb = n;

  blas_size_t info = 0;

  blaslapack::getrs(&trans, &n, &n_rhs, LILA_BLAS_CONST_CAST(coeff_t, A.data()),
                    &lda, LILA_BLAS_CONST_CAST(blas_size_t, ipiv.data()),
                    LILA_BLAS_CAST(coeff_t, B.data()), &ldb, &info);
  assert(info == 0);
}

template <class coeff_t>
inline void LUSolve(const Matrix<coeff_t> &A,
                    const std::vector<blas_size_t> &ipiv, Vector<coeff_t> &X,
                    char trans = 'N') {
  // check / get dimensions
  assert(A.nrows() == X.n());
  blas_size_t n = A.nrows();
  blas_size_t n_rhs = 1;
  blas_size_t lda = n;
  blas_size_t ldb = n;

  blas_size_t info = 0;

  blaslapack::getrs(&trans, &n, &n_rhs, LILA_BLAS_CONST_CAST(coeff_t, A.data()),
                    &lda, LILA_BLAS_CONST_CAST(blas_size_t, ipiv.data()),
                    LILA_BLAS_CAST(coeff_t, X.data()), &ldb, &info);
  assert(info == 0);
}

template <class coeff_t> inline void Invert(Matrix<coeff_t> &A) {

  // check / get dimensions
  assert(A.nrows() == A.ncols());
  blas_size_t n = A.nrows();

  std::vector<blas_size_t> ipiv(n);
  blas_size_t lwork = n * n;
  std::vector<coeff_t> work(lwork);
  blas_size_t info = 0;

  blaslapack::getrf(&n, &n, LILA_BLAS_CAST(coeff_t, A.data()), &n, ipiv.data(),
                    &info);
  assert(info == 0);
  blaslapack::getri(&n, LILA_BLAS_CAST(coeff_t, A.data()), &n, ipiv.data(),
                    LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);
  assert(info == 0);
}

template <class coeff_t>
inline std::vector<coeff_t> QRDecompose(Matrix<coeff_t> &A) {
  // check / get dimensions
  blas_size_t m = A.nrows();
  blas_size_t n = A.ncols();
  blas_size_t lda = m;

  std::vector<coeff_t> tau(std::min(m, n));
  blas_size_t info = 0;

  // get optimal work size
  blas_size_t lwork = -1;
  std::vector<coeff_t> work(1);
  blaslapack::geqrf(&m, &n, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                    LILA_BLAS_CAST(coeff_t, tau.data()),
                    LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);
  assert(info == 0);
  lwork = static_cast<blas_size_t>(real(work[0]));
  work.resize(lwork);

  // Run QR Decomposition
  blaslapack::geqrf(&m, &n, LILA_BLAS_CAST(coeff_t, A.data()), &lda,
                    LILA_BLAS_CAST(coeff_t, tau.data()),
                    LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);
  assert(info == 0);

  return tau;
}

template <class coeff_t>
inline Matrix<coeff_t> QRGetQ(Matrix<coeff_t> &A,
                              const std::vector<coeff_t> &tau) {
  // check / get dimensions
  blas_size_t m = A.nrows();
  blas_size_t k = tau.size();
  blas_size_t lda = m;

  Matrix<coeff_t> Q(A);
  blas_size_t info = 0;

  // get optimal work size
  blas_size_t lwork = -1;
  std::vector<coeff_t> work(1);
  blaslapack::orgqr(&m, &k, &k, LILA_BLAS_CAST(coeff_t, Q.data()), &lda,
                    LILA_BLAS_CONST_CAST(coeff_t, tau.data()),
                    LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);
  assert(info == 0);
  lwork = static_cast<blas_size_t>(real(work[0]));
  work.resize(lwork);

  // Compute Q
  blaslapack::orgqr(&m, &k, &k, LILA_BLAS_CAST(coeff_t, Q.data()), &lda,
                    LILA_BLAS_CONST_CAST(coeff_t, tau.data()),
                    LILA_BLAS_CAST(coeff_t, work.data()), &lwork, &info);

  Q.resize(m, k);
  assert(info == 0);
  return Q;
}

template <class coeff_t> inline Matrix<coeff_t> GetUpper(Matrix<coeff_t> &A) {
  blas_size_t m = A.nrows();
  blas_size_t n = A.ncols();
  blas_size_t k = std::min(m, n);

  Matrix<coeff_t> R(k, n);
  Zeros(R);
  for (int row = 0; row < k; ++row)
    for (int column = 0; column < n; ++column)
      if (row <= column)
        R(row, column) = A(row, column);
  return R;
}

} // namespace lila
