#pragma once

#include <complex>
#include <cstdlib>

#include "blaslapack_extern.h"
#include "blaslapack_types.h"

namespace lila::blaslapack {

// Copy
inline void copy(const blas_size_t *N, const blas_float_t *x,
                 const blas_size_t *incx, blas_float_t *y,
                 const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(scopy)(N, x, incx, y, incy);
}
inline void copy(const blas_size_t *N, const blas_double_t *x,
                 const blas_size_t *incx, blas_double_t *y,
                 const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(dcopy)(N, x, incx, y, incy);
}
inline void copy(const blas_size_t *N, const blas_scomplex_t *x,
                 const blas_size_t *incx, blas_scomplex_t *y,
                 const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(ccopy)(N, x, incx, y, incy);
}
inline void copy(const blas_size_t *N, const blas_complex_t *x,
                 const blas_size_t *incx, blas_complex_t *y,
                 const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(zcopy)(N, x, incx, y, incy);
}

// Axpy
inline void axpy(const blas_size_t *N, const blas_float_t *alpha,
                 const blas_float_t *x, const blas_size_t *incx,
                 blas_float_t *y, const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(saxpy)(N, alpha, x, incx, y, incy);
}
inline void axpy(const blas_size_t *N, const blas_double_t *alpha,
                 const blas_double_t *x, const blas_size_t *incx,
                 blas_double_t *y, const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(daxpy)(N, alpha, x, incx, y, incy);
}
inline void axpy(const blas_size_t *N, const blas_scomplex_t *alpha,
                 const blas_scomplex_t *x, const blas_size_t *incx,
                 blas_scomplex_t *y, const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(caxpy)(N, alpha, x, incx, y, incy);
}
inline void axpy(const blas_size_t *N, const blas_complex_t *alpha,
                 const blas_complex_t *x, const blas_size_t *incx,
                 blas_complex_t *y, const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(zaxpy)(N, alpha, x, incx, y, incy);
}

// Scal
inline void scal(const blas_size_t *N, const blas_float_t *alpha,
                 blas_float_t *x, const blas_size_t *incx) {
  __LAPACK_ROUTINE_NAME(sscal)(N, alpha, x, incx);
}
inline void scal(const blas_size_t *N, const blas_double_t *alpha,
                 blas_double_t *x, const blas_size_t *incx) {
  __LAPACK_ROUTINE_NAME(dscal)(N, alpha, x, incx);
}
inline void scal(const blas_size_t *N, const blas_scomplex_t *alpha,
                 blas_scomplex_t *x, const blas_size_t *incx) {
  __LAPACK_ROUTINE_NAME(cscal)(N, alpha, x, incx);
}
inline void scal(const blas_size_t *N, const blas_complex_t *alpha,
                 blas_complex_t *x, const blas_size_t *incx) {
  __LAPACK_ROUTINE_NAME(zscal)(N, alpha, x, incx);
}

// Dot
inline float dot(const blas_size_t *N, const blas_float_t *x,
                 const blas_size_t *incx, const blas_float_t *y,
                 const blas_size_t *incy) {
  return __LAPACK_ROUTINE_NAME(sdot)(N, x, incx, y, incy);
}
inline double dot(const blas_size_t *N, const blas_double_t *x,
                  const blas_size_t *incx, const blas_double_t *y,
                  const blas_size_t *incy) {
  return __LAPACK_ROUTINE_NAME(ddot)(N, x, incx, y, incy);
}

inline blas_scomplex_t dot(const blas_size_t *N, const blas_scomplex_t *x,
                           const blas_size_t *incx, const blas_scomplex_t *y,
                           const blas_size_t *incy) {
#if defined(LILA_USE_MKL) || defined(LILA_USE_ACCELERATE)
  blas_scomplex_t res;
  __LAPACK_ROUTINE_NAME(cdotc)(&res, N, x, incx, y, incy);
  return res;
#else
  return __LAPACK_ROUTINE_NAME(cdotc)(N, x, incx, y, incy);
#endif
}

inline blas_complex_t dot(const blas_size_t *N, const blas_complex_t *x,
                          const blas_size_t *incx, const blas_complex_t *y,
                          const blas_size_t *incy) {
#if defined(LILA_USE_MKL) || defined(LILA_USE_ACCELERATE)
  blas_complex_t res;
  __LAPACK_ROUTINE_NAME(zdotc)(&res, N, x, incx, y, incy);
  return res;
#else
  return __LAPACK_ROUTINE_NAME(zdotc)(N, x, incx, y, incy);
#endif
}

// Nrm2
inline blas_float_t nrm2(const blas_size_t *N, const blas_float_t *x,
                         const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(snrm2)(N, x, incx);
}
inline blas_double_t nrm2(const blas_size_t *N, const blas_double_t *x,
                          const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(dnrm2)(N, x, incx);
}
inline blas_float_t nrm2(const blas_size_t *N, const blas_scomplex_t *x,
                         const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(scnrm2)(N, x, incx);
}
inline blas_double_t nrm2(const blas_size_t *N, const blas_complex_t *x,
                          const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(dznrm2)(N, x, incx);
}

// Asum
inline blas_float_t asum(const blas_size_t *N, const blas_float_t *x,
                         const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(sasum)(N, x, incx);
}
inline blas_double_t asum(const blas_size_t *N, const blas_double_t *x,
                          const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(dasum)(N, x, incx);
}
inline blas_float_t asum(const blas_size_t *N, const blas_scomplex_t *x,
                         const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(scasum)(N, x, incx);
}
inline blas_double_t asum(const blas_size_t *N, const blas_complex_t *x,
                          const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(dzasum)(N, x, incx);
}

// IAmax
inline blas_size_t iamax(const blas_size_t *N, const blas_float_t *x,
                         const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(isamax)(N, x, incx);
}
inline blas_size_t iamax(const blas_size_t *N, const blas_double_t *x,
                         const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(idamax)(N, x, incx);
}
inline blas_size_t iamax(const blas_size_t *N, const blas_scomplex_t *x,
                         const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(icamax)(N, x, incx);
}
inline blas_size_t iamax(const blas_size_t *N, const blas_complex_t *x,
                         const blas_size_t *incx) {
  return __LAPACK_ROUTINE_NAME(izamax)(N, x, incx);
}

// Gemv
inline void gemv(const char *trans, const blas_size_t *m, const blas_size_t *n,
                 const blas_float_t *alpha, const blas_float_t *A,
                 const blas_size_t *dima, const blas_float_t *x,
                 const blas_size_t *incx, const blas_float_t *beta,
                 blas_float_t *y, const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(sgemv)
  (trans, m, n, alpha, A, dima, x, incx, beta, y, incy);
}
inline void gemv(const char *trans, const blas_size_t *m, const blas_size_t *n,
                 const blas_double_t *alpha, const blas_double_t *A,
                 const blas_size_t *dima, const blas_double_t *x,
                 const blas_size_t *incx, const blas_double_t *beta,
                 blas_double_t *y, const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(dgemv)
  (trans, m, n, alpha, A, dima, x, incx, beta, y, incy);
}
inline void gemv(const char *trans, const blas_size_t *m, const blas_size_t *n,
                 const blas_scomplex_t *alpha, const blas_scomplex_t *A,
                 const blas_size_t *dima, const blas_scomplex_t *x,
                 const blas_size_t *incx, const blas_scomplex_t *beta,
                 blas_scomplex_t *y, const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(cgemv)
  (trans, m, n, alpha, A, dima, x, incx, beta, y, incy);
}
inline void gemv(const char *trans, const blas_size_t *m, const blas_size_t *n,
                 const blas_complex_t *alpha, const blas_complex_t *A,
                 const blas_size_t *dima, const blas_complex_t *x,
                 const blas_size_t *incx, const blas_complex_t *beta,
                 blas_complex_t *y, const blas_size_t *incy) {
  __LAPACK_ROUTINE_NAME(zgemv)
  (trans, m, n, alpha, A, dima, x, incx, beta, y, incy);
}

// Gemm
inline void gemm(const char *transa, const char *transb, const blas_size_t *m,
                 const blas_size_t *n, const blas_size_t *k,
                 const blas_float_t *alpha, const blas_float_t *A,
                 const blas_size_t *dima, const blas_float_t *B,
                 const blas_size_t *dimb, const blas_float_t *beta,
                 blas_float_t *C, const blas_size_t *dimc) {
  __LAPACK_ROUTINE_NAME(sgemm)
  (transa, transb, m, n, k, alpha, A, dima, B, dimb, beta, C, dimc);
}
inline void gemm(const char *transa, const char *transb, const blas_size_t *m,
                 const blas_size_t *n, const blas_size_t *k,
                 const blas_double_t *alpha, const blas_double_t *A,
                 const blas_size_t *dima, const blas_double_t *B,
                 const blas_size_t *dimb, const blas_double_t *beta,
                 blas_double_t *C, const blas_size_t *dimc) {
  __LAPACK_ROUTINE_NAME(dgemm)
  (transa, transb, m, n, k, alpha, A, dima, B, dimb, beta, C, dimc);
}
inline void gemm(const char *transa, const char *transb, const blas_size_t *m,
                 const blas_size_t *n, const blas_size_t *k,
                 const blas_scomplex_t *alpha, const blas_scomplex_t *A,
                 const blas_size_t *dima, const blas_scomplex_t *B,
                 const blas_size_t *dimb, const blas_scomplex_t *beta,
                 blas_scomplex_t *C, const blas_size_t *dimc) {
  __LAPACK_ROUTINE_NAME(cgemm)
  (transa, transb, m, n, k, alpha, A, dima, B, dimb, beta, C, dimc);
}
inline void gemm(const char *transa, const char *transb, const blas_size_t *m,
                 const blas_size_t *n, const blas_size_t *k,
                 const blas_complex_t *alpha, const blas_complex_t *A,
                 const blas_size_t *dima, const blas_complex_t *B,
                 const blas_size_t *dimb, const blas_complex_t *beta,
                 blas_complex_t *C, const blas_size_t *dimc) {
  __LAPACK_ROUTINE_NAME(zgemm)
  (transa, transb, m, n, k, alpha, A, dima, B, dimb, beta, C, dimc);
}

//////////////////////////
// Linear Solve
// Gesv
inline void gesv(LAPACK_CONST blas_size_t *n, LAPACK_CONST blas_size_t *n_rhs,
                 blas_float_t *A, LAPACK_CONST blas_size_t *lda,
                 blas_size_t *ipiv, blas_float_t *B,
                 LAPACK_CONST blas_size_t *ldb, int *info) {
  __LAPACK_ROUTINE_NAME(sgesv)(n, n_rhs, A, lda, ipiv, B, ldb, info);
}
inline void gesv(LAPACK_CONST blas_size_t *n, LAPACK_CONST blas_size_t *n_rhs,
                 blas_double_t *A, LAPACK_CONST blas_size_t *lda,
                 blas_size_t *ipiv, blas_double_t *B,
                 LAPACK_CONST blas_size_t *ldb, int *info) {
  __LAPACK_ROUTINE_NAME(dgesv)(n, n_rhs, A, lda, ipiv, B, ldb, info);
}
inline void gesv(LAPACK_CONST blas_size_t *n, LAPACK_CONST blas_size_t *n_rhs,
                 blas_scomplex_t *A, LAPACK_CONST blas_size_t *lda,
                 blas_size_t *ipiv, blas_scomplex_t *B,
                 LAPACK_CONST blas_size_t *ldb, int *info) {
  __LAPACK_ROUTINE_NAME(cgesv)(n, n_rhs, A, lda, ipiv, B, ldb, info);
}
inline void gesv(LAPACK_CONST blas_size_t *n, LAPACK_CONST blas_size_t *n_rhs,
                 blas_complex_t *A, LAPACK_CONST blas_size_t *lda,
                 blas_size_t *ipiv, blas_complex_t *B,
                 LAPACK_CONST blas_size_t *ldb, int *info) {
  __LAPACK_ROUTINE_NAME(zgesv)(n, n_rhs, A, lda, ipiv, B, ldb, info);
}

//////////////////////////
// LU Decomposition
// Getrf (performs LU Decomposition)
inline void getrf(LAPACK_CONST blas_size_t *M, LAPACK_CONST blas_size_t *N,
                  blas_float_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_size_t *ipiv, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(sgetrf)(M, N, A, lda, ipiv, info);
}
inline void getrf(LAPACK_CONST blas_size_t *M, LAPACK_CONST blas_size_t *N,
                  blas_double_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_size_t *ipiv, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(dgetrf)(M, N, A, lda, ipiv, info);
}
inline void getrf(LAPACK_CONST blas_size_t *M, LAPACK_CONST blas_size_t *N,
                  blas_scomplex_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_size_t *ipiv, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(cgetrf)(M, N, A, lda, ipiv, info);
}
inline void getrf(LAPACK_CONST blas_size_t *M, LAPACK_CONST blas_size_t *N,
                  blas_complex_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_size_t *ipiv, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(zgetrf)(M, N, A, lda, ipiv, info);
}

// Getrs (solves system of equations from LU Decomposition)
inline void getrs(LAPACK_CONST char *trans, LAPACK_CONST blas_size_t *n,
                  LAPACK_CONST blas_size_t *n_rhs, LAPACK_CONST blas_float_t *A,
                  LAPACK_CONST blas_size_t *lda, LAPACK_CONST blas_size_t *ipiv,
                  blas_float_t *B, LAPACK_CONST blas_size_t *ldb,
                  blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(sgetrs)(trans, n, n_rhs, A, lda, ipiv, B, ldb, info);
}
inline void getrs(LAPACK_CONST char *trans, LAPACK_CONST blas_size_t *n,
                  LAPACK_CONST blas_size_t *n_rhs,
                  LAPACK_CONST blas_double_t *A, LAPACK_CONST blas_size_t *lda,
                  LAPACK_CONST blas_size_t *ipiv, blas_double_t *B,
                  LAPACK_CONST blas_size_t *ldb, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(dgetrs)(trans, n, n_rhs, A, lda, ipiv, B, ldb, info);
}
inline void getrs(LAPACK_CONST char *trans, LAPACK_CONST blas_size_t *n,
                  LAPACK_CONST blas_size_t *n_rhs,
                  LAPACK_CONST blas_scomplex_t *A,
                  LAPACK_CONST blas_size_t *lda, LAPACK_CONST blas_size_t *ipiv,
                  blas_scomplex_t *B, LAPACK_CONST blas_size_t *ldb,
                  blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(cgetrs)(trans, n, n_rhs, A, lda, ipiv, B, ldb, info);
}
inline void getrs(LAPACK_CONST char *trans, LAPACK_CONST blas_size_t *n,
                  LAPACK_CONST blas_size_t *n_rhs,
                  LAPACK_CONST blas_complex_t *A, LAPACK_CONST blas_size_t *lda,
                  LAPACK_CONST blas_size_t *ipiv, blas_complex_t *B,
                  LAPACK_CONST blas_size_t *ldb, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(zgetrs)(trans, n, n_rhs, A, lda, ipiv, B, ldb, info);
}

// Getri (computes inverse from LU Decomposition)
inline void getri(LAPACK_CONST blas_size_t *N, blas_float_t *A,
                  LAPACK_CONST blas_size_t *lda, LAPACK_CONST blas_size_t *IPIV,
                  float *WORK, LAPACK_CONST blas_size_t *lwork,
                  blas_size_t *INFO) {
  __LAPACK_ROUTINE_NAME(sgetri)(N, A, lda, IPIV, WORK, lwork, INFO);
}
inline void getri(LAPACK_CONST blas_size_t *N, blas_double_t *A,
                  LAPACK_CONST blas_size_t *lda, LAPACK_CONST blas_size_t *IPIV,
                  blas_double_t *WORK, LAPACK_CONST blas_size_t *lwork,
                  blas_size_t *INFO) {
  __LAPACK_ROUTINE_NAME(dgetri)(N, A, lda, IPIV, WORK, lwork, INFO);
}
inline void getri(LAPACK_CONST blas_size_t *N, blas_scomplex_t *A,
                  LAPACK_CONST blas_size_t *lda, LAPACK_CONST blas_size_t *IPIV,
                  blas_scomplex_t *WORK, LAPACK_CONST blas_size_t *lwork,
                  blas_size_t *INFO) {
  __LAPACK_ROUTINE_NAME(cgetri)(N, A, lda, IPIV, WORK, lwork, INFO);
}
inline void getri(LAPACK_CONST blas_size_t *N, blas_complex_t *A,
                  LAPACK_CONST blas_size_t *lda, LAPACK_CONST blas_size_t *IPIV,
                  blas_complex_t *WORK, LAPACK_CONST blas_size_t *lwork,
                  blas_size_t *INFO) {
  __LAPACK_ROUTINE_NAME(zgetri)(N, A, lda, IPIV, WORK, lwork, INFO);
}

//////////////////////////
// QR Decomposition
// Geqrf  (does the QR Decomposition)
inline void geqrf(LAPACK_CONST blas_size_t *m, LAPACK_CONST blas_size_t *n,
                  blas_float_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_float_t *tau, blas_float_t *work,
                  LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(sgeqrf)(m, n, A, lda, tau, work, lwork, info);
}
inline void geqrf(LAPACK_CONST blas_size_t *m, LAPACK_CONST blas_size_t *n,
                  blas_double_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_double_t *tau, blas_double_t *work,
                  LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(dgeqrf)(m, n, A, lda, tau, work, lwork, info);
}
inline void geqrf(LAPACK_CONST blas_size_t *m, LAPACK_CONST blas_size_t *n,
                  blas_scomplex_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_scomplex_t *tau, blas_scomplex_t *work,
                  LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(cgeqrf)(m, n, A, lda, tau, work, lwork, info);
}
inline void geqrf(LAPACK_CONST blas_size_t *m, LAPACK_CONST blas_size_t *n,
                  blas_complex_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_complex_t *tau, blas_complex_t *work,
                  LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(zgeqrf)(m, n, A, lda, tau, work, lwork, info);
}

// Orgqr / Ungqr (retrieves the matrix Q of the QR Decomposition)
// (Note: use orgqr also for the complex routines xungqr)
inline void orgqr(LAPACK_CONST blas_size_t *m, LAPACK_CONST blas_size_t *n,
                  LAPACK_CONST blas_size_t *k, blas_float_t *A,
                  LAPACK_CONST blas_size_t *lda, LAPACK_CONST blas_float_t *tau,
                  blas_float_t *work, LAPACK_CONST blas_size_t *lwork,
                  blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(sorgqr)(m, n, k, A, lda, tau, work, lwork, info);
}
inline void orgqr(LAPACK_CONST blas_size_t *m, LAPACK_CONST blas_size_t *n,
                  LAPACK_CONST blas_size_t *k, blas_double_t *A,
                  LAPACK_CONST blas_size_t *lda,
                  LAPACK_CONST blas_double_t *tau, blas_double_t *work,
                  LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(dorgqr)(m, n, k, A, lda, tau, work, lwork, info);
}
inline void orgqr(LAPACK_CONST blas_size_t *m, LAPACK_CONST blas_size_t *n,
                  LAPACK_CONST blas_size_t *k, blas_scomplex_t *A,
                  LAPACK_CONST blas_size_t *lda,
                  LAPACK_CONST blas_scomplex_t *tau, blas_scomplex_t *work,
                  LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(cungqr)(m, n, k, A, lda, tau, work, lwork, info);
}
inline void orgqr(LAPACK_CONST blas_size_t *m, LAPACK_CONST blas_size_t *n,
                  LAPACK_CONST blas_size_t *k, blas_complex_t *A,
                  LAPACK_CONST blas_size_t *lda,
                  LAPACK_CONST blas_complex_t *tau, blas_complex_t *work,
                  LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(zungqr)(m, n, k, A, lda, tau, work, lwork, info);
}

//////////////////////////
// Cholesky Decomposition
inline void potrf(LAPACK_CONST char *uplo, LAPACK_CONST blas_size_t *n,
                  blas_float_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(spotrf)(uplo, n, A, lda, info);
}

inline void potrf(LAPACK_CONST char *uplo, LAPACK_CONST blas_size_t *n,
                  blas_double_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(dpotrf)(uplo, n, A, lda, info);
}

inline void potrf(LAPACK_CONST char *uplo, LAPACK_CONST blas_size_t *n,
                  blas_scomplex_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(cpotrf)(uplo, n, A, lda, info);
}

inline void potrf(LAPACK_CONST char *uplo, LAPACK_CONST blas_size_t *n,
                  blas_complex_t *A, LAPACK_CONST blas_size_t *lda,
                  blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(zpotrf)(uplo, n, A, lda, info);
}

//////////////////////////
// Eigenvalues

// symmetric/hermitian
inline void syev(LAPACK_CONST char *jobz, LAPACK_CONST char *uplo,
                 LAPACK_CONST blas_size_t *n, blas_float_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_float_t *w,
                 blas_float_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(ssyev)(jobz, uplo, n, a, lda, w, work, lwork, info);
}

inline void syev(LAPACK_CONST char *jobz, LAPACK_CONST char *uplo,
                 LAPACK_CONST blas_size_t *n, blas_double_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_double_t *w,
                 blas_double_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(dsyev)(jobz, uplo, n, a, lda, w, work, lwork, info);
}

inline void heev(LAPACK_CONST char *jobz, LAPACK_CONST char *uplo,
                 LAPACK_CONST blas_size_t *n, blas_scomplex_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_float_t *w,
                 blas_scomplex_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_float_t *rwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(cheev)
  (jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

inline void heev(LAPACK_CONST char *jobz, LAPACK_CONST char *uplo,
                 LAPACK_CONST blas_size_t *n, blas_complex_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_double_t *w,
                 blas_complex_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_double_t *rwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(zheev)
  (jobz, uplo, n, a, lda, w, work, lwork, rwork, info);
}

// (not actual lapack routines, create extra rwork)
inline void syev(LAPACK_CONST char *jobz, LAPACK_CONST char *uplo,
                 LAPACK_CONST blas_size_t *n, blas_scomplex_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_float_t *w,
                 blas_scomplex_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_size_t *info) {
  std::vector<blas_float_t> rwork(3 * (*n) - 2);
  heev(jobz, uplo, n, a, lda, w, work, lwork, rwork.data(), info);
}

inline void syev(LAPACK_CONST char *jobz, LAPACK_CONST char *uplo,
                 LAPACK_CONST blas_size_t *n, blas_complex_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_double_t *w,
                 blas_complex_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_size_t *info) {
  std::vector<blas_double_t> rwork(3 * (*n) - 2);
  heev(jobz, uplo, n, a, lda, w, work, lwork, rwork.data(), info);
}

// generic
inline void geev(LAPACK_CONST char *jobvl, LAPACK_CONST char *jobvr,
                 LAPACK_CONST blas_size_t *n, blas_float_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_float_t *wr,
                 blas_float_t *wi, blas_float_t *vl,
                 LAPACK_CONST blas_size_t *ldvl, blas_float_t *vr,
                 LAPACK_CONST blas_size_t *ldvr, blas_float_t *work,
                 LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(sgeev)
  (jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}
inline void geev(LAPACK_CONST char *jobvl, LAPACK_CONST char *jobvr,
                 LAPACK_CONST blas_size_t *n, blas_double_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_double_t *wr,
                 blas_double_t *wi, blas_double_t *vl,
                 LAPACK_CONST blas_size_t *ldvl, blas_double_t *vr,
                 LAPACK_CONST blas_size_t *ldvr, blas_double_t *work,
                 LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(dgeev)
  (jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
}

// (not actual lapack routines, convert real/imag vectors to complex)
inline void geev(LAPACK_CONST char *jobvl, LAPACK_CONST char *jobvr,
                 LAPACK_CONST blas_size_t *n, blas_float_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_scomplex_t *w,
                 blas_float_t *vl, LAPACK_CONST blas_size_t *ldvl,
                 blas_float_t *vr, LAPACK_CONST blas_size_t *ldvr,
                 blas_float_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_size_t *info) {
  std::vector<blas_float_t> wr(*n);
  std::vector<blas_float_t> wi(*n);
  geev(jobvl, jobvr, n, a, lda, wr.data(), wi.data(), vl, ldvl, vr, ldvr, work,
       lwork, info);
  for (blas_size_t i = 0; i < *n; ++i) {
    // w[i] = blas_scomplex_t(wr[i], wi[i]);
    w[i] = {wr[i], wi[i]};
  }
}

inline void geev(LAPACK_CONST char *jobvl, LAPACK_CONST char *jobvr,
                 LAPACK_CONST blas_size_t *n, blas_double_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_complex_t *w,
                 blas_double_t *vl, LAPACK_CONST blas_size_t *ldvl,
                 blas_double_t *vr, LAPACK_CONST blas_size_t *ldvr,
                 blas_double_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_size_t *info) {
  std::vector<blas_double_t> wr(*n);
  std::vector<blas_double_t> wi(*n);
  geev(jobvl, jobvr, n, a, lda, wr.data(), wi.data(), vl, ldvl, vr, ldvr, work,
       lwork, info);
  for (blas_size_t i = 0; i < *n; ++i) {
    // w[i] = blas_complex_t(wr[i], wi[i]);
    w[i] = {wr[i], wi[i]};
  }
}

inline void geev(LAPACK_CONST char *jobvl, LAPACK_CONST char *jobvr,
                 LAPACK_CONST blas_size_t *n, blas_scomplex_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_scomplex_t *w,
                 blas_scomplex_t *vl, LAPACK_CONST blas_size_t *ldvl,
                 blas_scomplex_t *vr, LAPACK_CONST blas_size_t *ldvr,
                 blas_scomplex_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_size_t *info) {
  std::vector<blas_float_t> rwork(std::max(1, *n * 2));
  __LAPACK_ROUTINE_NAME(cgeev)
  (jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork.data(),
   info);
}

inline void geev(LAPACK_CONST char *jobvl, LAPACK_CONST char *jobvr,
                 LAPACK_CONST blas_size_t *n, blas_complex_t *a,
                 LAPACK_CONST blas_size_t *lda, blas_complex_t *w,
                 blas_complex_t *vl, LAPACK_CONST blas_size_t *ldvl,
                 blas_complex_t *vr, LAPACK_CONST blas_size_t *ldvr,
                 blas_complex_t *work, LAPACK_CONST blas_size_t *lwork,
                 blas_size_t *info) {
  std::vector<blas_double_t> rwork(std::max(1, *n * 2));
  __LAPACK_ROUTINE_NAME(zgeev)
  (jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork.data(),
   info);
}

// Real symmetric tridiagonal eigensolvers
inline void stev(LAPACK_CONST char *jobz, LAPACK_CONST blas_size_t *N,
                 blas_float_t *D, blas_float_t *E, blas_float_t *Z,
                 LAPACK_CONST blas_size_t *ldz, blas_float_t *work,
                 blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(sstev)(jobz, N, D, E, Z, ldz, work, info);
}

inline void stev(LAPACK_CONST char *jobz, LAPACK_CONST blas_size_t *N,
                 blas_double_t *D, blas_double_t *E, blas_double_t *Z,
                 LAPACK_CONST blas_size_t *ldz, blas_double_t *work,
                 blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(dstev)(jobz, N, D, E, Z, ldz, work, info);
}

/////////////////////////////////////////////////////////////////////////
// Generalized eigenvalue problems
// (use sy.. instead of he.. for complex routines)
inline void sygv(LAPACK_CONST blas_size_t *itype, LAPACK_CONST char *jobz,
                 LAPACK_CONST char *uplo, LAPACK_CONST blas_size_t *N,
                 blas_float_t *A, LAPACK_CONST blas_size_t *lda,
                 blas_float_t *B, LAPACK_CONST blas_size_t *ldb,
                 blas_float_t *W, blas_float_t *work,
                 LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(ssygv)
  (itype, jobz, uplo, N, A, lda, B, ldb, W, work, lwork, info);
}

inline void sygv(LAPACK_CONST blas_size_t *itype, LAPACK_CONST char *jobz,
                 LAPACK_CONST char *uplo, LAPACK_CONST blas_size_t *N,
                 blas_double_t *A, LAPACK_CONST blas_size_t *lda,
                 blas_double_t *B, LAPACK_CONST blas_size_t *ldb,
                 blas_double_t *W, blas_double_t *work,
                 LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  __LAPACK_ROUTINE_NAME(dsygv)
  (itype, jobz, uplo, N, A, lda, B, ldb, W, work, lwork, info);
}

inline void sygv(LAPACK_CONST blas_size_t *itype, LAPACK_CONST char *jobz,
                 LAPACK_CONST char *uplo, LAPACK_CONST blas_size_t *N,
                 blas_scomplex_t *A, LAPACK_CONST blas_size_t *lda,
                 blas_scomplex_t *B, LAPACK_CONST blas_size_t *ldb,
                 blas_float_t *W, blas_scomplex_t *work,
                 LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  std::vector<blas_float_t> rwork(std::max(1, *N * 3 - 2));
  __LAPACK_ROUTINE_NAME(chegv)
  (itype, jobz, uplo, N, A, lda, B, ldb, W, work, lwork, rwork.data(), info);
}

inline void sygv(LAPACK_CONST blas_size_t *itype, LAPACK_CONST char *jobz,
                 LAPACK_CONST char *uplo, LAPACK_CONST blas_size_t *N,
                 blas_complex_t *A, LAPACK_CONST blas_size_t *lda,
                 blas_complex_t *B, LAPACK_CONST blas_size_t *ldb,
                 blas_double_t *W, blas_complex_t *work,
                 LAPACK_CONST blas_size_t *lwork, blas_size_t *info) {
  std::vector<blas_double_t> rwork(std::max(1, *N * 3 - 2));
  __LAPACK_ROUTINE_NAME(zhegv)
  (itype, jobz, uplo, N, A, lda, B, ldb, W, work, lwork, rwork.data(), info);
}

} // namespace lila::blaslapack
