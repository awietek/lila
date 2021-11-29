#pragma once

#include "blaslapack_types.h"

using blas_size_t = lila::blas_size_t;
using blas_float_t = lila::blas_float_t;
using blas_double_t = lila::blas_double_t;
using blas_scomplex_t = lila::blas_scomplex_t;
using blas_complex_t = lila::blas_complex_t;
using lapack_ret_t = lila::lapack_ret_t;

#if defined(LILA_USE_LAPACK) // or defined(LILA_USE_ACCELERATE)

// Copy
extern "C" void scopy_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       blas_float_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void dcopy_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       blas_double_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void ccopy_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       blas_scomplex_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void zcopy_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       blas_complex_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);

// Axpy
extern "C" void saxpy_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       blas_float_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void daxpy_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       blas_double_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void caxpy_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       blas_scomplex_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void zaxpy_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       blas_complex_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);

// Scal
extern "C" void sscal_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *alpha,
                       blas_float_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" void dscal_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *alpha,
                       blas_double_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" void cscal_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *alpha,
                       blas_scomplex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" void zscal_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *alpha,
                       blas_complex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx);

// Dot
extern "C" float sdot_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" double ddot_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                        __LILA_BLAS_LAPACK_CONST blas_double_t *x,
                        __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                        __LILA_BLAS_LAPACK_CONST blas_double_t *y,
                        __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
#ifdef LILA_USE_ACCELERATE
extern "C" void cdotc_(blas_scomplex_t *ret,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void zdotc_(blas_complex_t *ret,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *N,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
#else
extern "C" blas_scomplex_t cdotc_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                  __LILA_BLAS_LAPACK_CONST blas_scomplex_t *x,
                                  __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                                  __LILA_BLAS_LAPACK_CONST blas_scomplex_t *y,
                                  __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" blas_complex_t zdotc_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                 __LILA_BLAS_LAPACK_CONST blas_complex_t *x,
                                 __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                                 __LILA_BLAS_LAPACK_CONST blas_complex_t *y,
                                 __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
#endif

// Norm2
extern "C" blas_float_t snrm2_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                               __LILA_BLAS_LAPACK_CONST blas_float_t *x,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" blas_double_t dnrm2_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                __LILA_BLAS_LAPACK_CONST blas_double_t *x,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" blas_float_t scnrm2_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                __LILA_BLAS_LAPACK_CONST blas_scomplex_t *x,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" blas_double_t dznrm2_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                 __LILA_BLAS_LAPACK_CONST blas_complex_t *x,
                                 __LILA_BLAS_LAPACK_CONST blas_size_t *incx);

// Asum
extern "C" blas_float_t sasum_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                               __LILA_BLAS_LAPACK_CONST blas_float_t *x,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" blas_double_t dasum_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                __LILA_BLAS_LAPACK_CONST blas_double_t *x,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" blas_float_t scasum_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                __LILA_BLAS_LAPACK_CONST blas_scomplex_t *x,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" blas_double_t dzasum_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                 __LILA_BLAS_LAPACK_CONST blas_complex_t *x,
                                 __LILA_BLAS_LAPACK_CONST blas_size_t *incx);

// IAmax
extern "C" blas_size_t isamax_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                               __LILA_BLAS_LAPACK_CONST blas_float_t *x,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" blas_size_t idamax_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                               __LILA_BLAS_LAPACK_CONST blas_double_t *x,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" blas_size_t icamax_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                               __LILA_BLAS_LAPACK_CONST blas_scomplex_t *x,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *incx);
extern "C" blas_size_t izamax_(__LILA_BLAS_LAPACK_CONST blas_size_t *N,
                               __LILA_BLAS_LAPACK_CONST blas_complex_t *x,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *incx);

// Gemv
extern "C" void sgemv_(__LILA_BLAS_LAPACK_CONST char *trans,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *m,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *A,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dima,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *beta,
                       blas_float_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void dgemv_(__LILA_BLAS_LAPACK_CONST char *trans,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *m,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *A,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dima,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *beta,
                       blas_double_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void cgemv_(__LILA_BLAS_LAPACK_CONST char *trans,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *m,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *A,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dima,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *beta,
                       blas_scomplex_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);
extern "C" void zgemv_(__LILA_BLAS_LAPACK_CONST char *trans,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *m,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *A,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dima,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *x,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incx,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *beta,
                       blas_complex_t *y,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *incy);

// Gemm
extern "C" void sgemm_(__LILA_BLAS_LAPACK_CONST char *transa,
                       __LILA_BLAS_LAPACK_CONST char *transb,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *m,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *k,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *A,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dima,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *B,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dimb,
                       __LILA_BLAS_LAPACK_CONST blas_float_t *beta,
                       blas_float_t *C,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dimc);
extern "C" void dgemm_(__LILA_BLAS_LAPACK_CONST char *transa,
                       __LILA_BLAS_LAPACK_CONST char *transb,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *m,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *k,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *A,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dima,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *B,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dimb,
                       __LILA_BLAS_LAPACK_CONST blas_double_t *beta,
                       blas_double_t *C,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dimc);
extern "C" void cgemm_(__LILA_BLAS_LAPACK_CONST char *transa,
                       __LILA_BLAS_LAPACK_CONST char *transb,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *m,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *k,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *A,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dima,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *B,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dimb,
                       __LILA_BLAS_LAPACK_CONST blas_scomplex_t *beta,
                       blas_scomplex_t *C,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dimc);
extern "C" void zgemm_(__LILA_BLAS_LAPACK_CONST char *transa,
                       __LILA_BLAS_LAPACK_CONST char *transb,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *m,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *k,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *alpha,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *A,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dima,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *B,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dimb,
                       __LILA_BLAS_LAPACK_CONST blas_complex_t *beta,
                       blas_complex_t *C,
                       __LILA_BLAS_LAPACK_CONST blas_size_t *dimc);
#endif

#ifdef LILA_USE_LAPACK
//////////////////////////
// Linear Solve
// Gesv
extern "C" lapack_ret_t
sgesv_(__LILA_BLAS_LAPACK_CONST blas_size_t *n,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n_rhs, blas_float_t *A,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_size_t *ipiv,
       blas_float_t *B, __LILA_BLAS_LAPACK_CONST blas_size_t *ldb, int *info);
extern "C" lapack_ret_t
dgesv_(__LILA_BLAS_LAPACK_CONST blas_size_t *n,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n_rhs, blas_double_t *A,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_size_t *ipiv,
       blas_double_t *B, __LILA_BLAS_LAPACK_CONST blas_size_t *ldb, int *info);
extern "C" lapack_ret_t cgesv_(__LILA_BLAS_LAPACK_CONST blas_size_t *n,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *n_rhs,
                               blas_scomplex_t *A,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                               blas_size_t *ipiv, blas_scomplex_t *B,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *ldb,
                               int *info);
extern "C" lapack_ret_t
zgesv_(__LILA_BLAS_LAPACK_CONST blas_size_t *n,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n_rhs, blas_complex_t *A,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_size_t *ipiv,
       blas_complex_t *B, __LILA_BLAS_LAPACK_CONST blas_size_t *ldb, int *info);

//////////////////////////
// LU Decomposition
// Getrf (performs LU Decomposition)
extern "C" lapack_ret_t sgetrf_(__LILA_BLAS_LAPACK_CONST blas_size_t *M,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                blas_float_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_size_t *ipiv, blas_size_t *info);
extern "C" lapack_ret_t dgetrf_(__LILA_BLAS_LAPACK_CONST blas_size_t *M,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                blas_double_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_size_t *ipiv, blas_size_t *info);
extern "C" lapack_ret_t cgetrf_(__LILA_BLAS_LAPACK_CONST blas_size_t *M,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                blas_scomplex_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_size_t *ipiv, blas_size_t *info);
extern "C" lapack_ret_t zgetrf_(__LILA_BLAS_LAPACK_CONST blas_size_t *M,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *N,
                                blas_complex_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_size_t *ipiv, blas_size_t *info);

// Getrs (solves system of equations from LU Decomposition)
extern "C" lapack_ret_t sgetrs_(__LILA_BLAS_LAPACK_CONST char *trans,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n_rhs,
                                __LILA_BLAS_LAPACK_CONST blas_float_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *ipiv,
                                blas_float_t *B,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *ldb,
                                blas_size_t *info);
extern "C" lapack_ret_t dgetrs_(__LILA_BLAS_LAPACK_CONST char *trans,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n_rhs,
                                __LILA_BLAS_LAPACK_CONST blas_double_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *ipiv,
                                blas_double_t *B,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *ldb,
                                blas_size_t *info);
extern "C" lapack_ret_t cgetrs_(__LILA_BLAS_LAPACK_CONST char *trans,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n_rhs,
                                __LILA_BLAS_LAPACK_CONST blas_scomplex_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *ipiv,
                                blas_scomplex_t *B,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *ldb,
                                blas_size_t *info);
extern "C" lapack_ret_t zgetrs_(__LILA_BLAS_LAPACK_CONST char *trans,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n_rhs,
                                __LILA_BLAS_LAPACK_CONST blas_complex_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *ipiv,
                                blas_complex_t *B,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *ldb,
                                blas_size_t *info);

// Getri (computes inverse from LU Decomposition)
extern "C" lapack_ret_t
sgetri_(__LILA_BLAS_LAPACK_CONST blas_size_t *N, blas_float_t *A,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
        __LILA_BLAS_LAPACK_CONST blas_size_t *IPIV, blas_float_t *WORK,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *INFO);
extern "C" lapack_ret_t
dgetri_(__LILA_BLAS_LAPACK_CONST blas_size_t *N, blas_double_t *A,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
        __LILA_BLAS_LAPACK_CONST blas_size_t *IPIV, blas_double_t *WORK,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *INFO);
extern "C" lapack_ret_t
cgetri_(__LILA_BLAS_LAPACK_CONST blas_size_t *N, blas_scomplex_t *A,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
        __LILA_BLAS_LAPACK_CONST blas_size_t *IPIV, blas_scomplex_t *WORK,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *INFO);
extern "C" lapack_ret_t
zgetri_(__LILA_BLAS_LAPACK_CONST blas_size_t *N, blas_complex_t *A,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
        __LILA_BLAS_LAPACK_CONST blas_size_t *IPIV, blas_complex_t *WORK,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *INFO);

//////////////////////////
// QR Decomposition
// Geqrf (does the QR Decomposition)
extern "C" lapack_ret_t sgeqrf_(__LILA_BLAS_LAPACK_CONST blas_size_t *m,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                blas_float_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_float_t *tau, blas_float_t *work,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
                                blas_size_t *info);
extern "C" lapack_ret_t dgeqrf_(__LILA_BLAS_LAPACK_CONST blas_size_t *m,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                blas_double_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_double_t *tau, blas_double_t *work,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
                                blas_size_t *info);
extern "C" lapack_ret_t cgeqrf_(__LILA_BLAS_LAPACK_CONST blas_size_t *m,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                blas_scomplex_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_scomplex_t *tau, blas_scomplex_t *work,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
                                blas_size_t *info);
extern "C" lapack_ret_t zgeqrf_(__LILA_BLAS_LAPACK_CONST blas_size_t *m,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                blas_complex_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_complex_t *tau, blas_complex_t *work,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
                                blas_size_t *info);

// Orgqr / Ungqr (retrieves the matrix Q of the QR Decomposition)
extern "C" lapack_ret_t
sorgqr_(__LILA_BLAS_LAPACK_CONST blas_size_t *m,
        __LILA_BLAS_LAPACK_CONST blas_size_t *n,
        __LILA_BLAS_LAPACK_CONST blas_size_t *k, blas_float_t *A,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
        __LILA_BLAS_LAPACK_CONST blas_float_t *tau, blas_float_t *work,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *info);
extern "C" lapack_ret_t
dorgqr_(__LILA_BLAS_LAPACK_CONST blas_size_t *m,
        __LILA_BLAS_LAPACK_CONST blas_size_t *n,
        __LILA_BLAS_LAPACK_CONST blas_size_t *k, blas_double_t *A,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
        __LILA_BLAS_LAPACK_CONST blas_double_t *tau, blas_double_t *work,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *info);
extern "C" lapack_ret_t
cungqr_(__LILA_BLAS_LAPACK_CONST blas_size_t *m,
        __LILA_BLAS_LAPACK_CONST blas_size_t *n,
        __LILA_BLAS_LAPACK_CONST blas_size_t *k, blas_scomplex_t *A,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
        __LILA_BLAS_LAPACK_CONST blas_scomplex_t *tau, blas_scomplex_t *work,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *info);
extern "C" lapack_ret_t
zungqr_(__LILA_BLAS_LAPACK_CONST blas_size_t *m,
        __LILA_BLAS_LAPACK_CONST blas_size_t *n,
        __LILA_BLAS_LAPACK_CONST blas_size_t *k, blas_complex_t *A,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
        __LILA_BLAS_LAPACK_CONST blas_complex_t *tau, blas_complex_t *work,
        __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *info);

//////////////////////////
// Cholesky Decomposition
extern "C" lapack_ret_t spotrf_(__LILA_BLAS_LAPACK_CONST char *uplo,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                blas_float_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_size_t *info);

extern "C" lapack_ret_t dpotrf_(__LILA_BLAS_LAPACK_CONST char *uplo,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                blas_double_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_size_t *info);

extern "C" lapack_ret_t cpotrf_(__LILA_BLAS_LAPACK_CONST char *uplo,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                blas_scomplex_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_size_t *info);

extern "C" lapack_ret_t zpotrf_(__LILA_BLAS_LAPACK_CONST char *uplo,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *n,
                                blas_complex_t *A,
                                __LILA_BLAS_LAPACK_CONST blas_size_t *lda,
                                blas_size_t *info);

//////////////////////////
// Eigenvalues

// symmetric/hermitian
extern "C" lapack_ret_t
ssyev_(__LILA_BLAS_LAPACK_CONST char *jobz, __LILA_BLAS_LAPACK_CONST char *uplo,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n, blas_float_t *a,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_float_t *w,
       blas_float_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_size_t *info);

extern "C" lapack_ret_t
dsyev_(__LILA_BLAS_LAPACK_CONST char *jobz, __LILA_BLAS_LAPACK_CONST char *uplo,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n, blas_double_t *a,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_double_t *w,
       blas_double_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_size_t *info);

extern "C" lapack_ret_t
cheev_(__LILA_BLAS_LAPACK_CONST char *jobz, __LILA_BLAS_LAPACK_CONST char *uplo,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n, blas_scomplex_t *a,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_float_t *w,
       blas_scomplex_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_float_t *rwork, blas_size_t *info);

extern "C" lapack_ret_t
zheev_(__LILA_BLAS_LAPACK_CONST char *jobz, __LILA_BLAS_LAPACK_CONST char *uplo,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n, blas_complex_t *a,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_double_t *w,
       blas_complex_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_double_t *rwork, blas_size_t *info);

// generic
extern "C" lapack_ret_t
sgeev_(__LILA_BLAS_LAPACK_CONST char *jobvl,
       __LILA_BLAS_LAPACK_CONST char *jobvr,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n, blas_float_t *a,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_float_t *wr,
       blas_float_t *wi, blas_float_t *vl,
       __LILA_BLAS_LAPACK_CONST blas_size_t *ldvl, blas_float_t *vr,
       __LILA_BLAS_LAPACK_CONST blas_size_t *ldvr, blas_float_t *work,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *info);
extern "C" lapack_ret_t
dgeev_(__LILA_BLAS_LAPACK_CONST char *jobvl,
       __LILA_BLAS_LAPACK_CONST char *jobvr,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n, blas_double_t *a,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_double_t *wr,
       blas_double_t *wi, blas_double_t *vl,
       __LILA_BLAS_LAPACK_CONST blas_size_t *ldvl, blas_double_t *vr,
       __LILA_BLAS_LAPACK_CONST blas_size_t *ldvr, blas_double_t *work,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lwork, blas_size_t *info);
extern "C" lapack_ret_t
cgeev_(__LILA_BLAS_LAPACK_CONST char *jobvl,
       __LILA_BLAS_LAPACK_CONST char *jobvr,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n, blas_scomplex_t *a,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_scomplex_t *w,
       blas_scomplex_t *vl, __LILA_BLAS_LAPACK_CONST blas_size_t *ldvl,
       blas_scomplex_t *vr, __LILA_BLAS_LAPACK_CONST blas_size_t *ldvr,
       blas_scomplex_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_float_t *rwork, blas_size_t *info);
extern "C" lapack_ret_t
zgeev_(__LILA_BLAS_LAPACK_CONST char *jobvl,
       __LILA_BLAS_LAPACK_CONST char *jobvr,
       __LILA_BLAS_LAPACK_CONST blas_size_t *n, blas_complex_t *a,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_complex_t *w,
       blas_complex_t *vl, __LILA_BLAS_LAPACK_CONST blas_size_t *ldvl,
       blas_complex_t *vr, __LILA_BLAS_LAPACK_CONST blas_size_t *ldvr,
       blas_complex_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_double_t *rwork, blas_size_t *info);

// Real symmetric tridiagonal eigensolvers
extern "C" lapack_ret_t sstev_(__LILA_BLAS_LAPACK_CONST char *jobz,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *N,
                               blas_float_t *D, blas_float_t *E,
                               blas_float_t *Z,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *ldz,
                               blas_float_t *work, blas_size_t *info);

extern "C" lapack_ret_t dstev_(__LILA_BLAS_LAPACK_CONST char *jobz,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *N,
                               blas_double_t *D, blas_double_t *E,
                               blas_double_t *Z,
                               __LILA_BLAS_LAPACK_CONST blas_size_t *ldz,
                               blas_double_t *work, blas_size_t *info);

///////////////////////////////////
// Generalized eigenvalue problems
extern "C" lapack_ret_t
ssygv_(__LILA_BLAS_LAPACK_CONST blas_size_t *itype,
       __LILA_BLAS_LAPACK_CONST char *jobz, __LILA_BLAS_LAPACK_CONST char *uplo,
       __LILA_BLAS_LAPACK_CONST blas_size_t *N, blas_float_t *A,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_float_t *B,
       __LILA_BLAS_LAPACK_CONST blas_size_t *ldb, blas_float_t *W,
       blas_float_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_size_t *info);

extern "C" lapack_ret_t
dsygv_(__LILA_BLAS_LAPACK_CONST blas_size_t *itype,
       __LILA_BLAS_LAPACK_CONST char *jobz, __LILA_BLAS_LAPACK_CONST char *uplo,
       __LILA_BLAS_LAPACK_CONST blas_size_t *N, blas_double_t *A,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_double_t *B,
       __LILA_BLAS_LAPACK_CONST blas_size_t *ldb, blas_double_t *W,
       blas_double_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_size_t *info);

extern "C" lapack_ret_t
chegv_(__LILA_BLAS_LAPACK_CONST blas_size_t *itype,
       __LILA_BLAS_LAPACK_CONST char *jobz, __LILA_BLAS_LAPACK_CONST char *uplo,
       __LILA_BLAS_LAPACK_CONST blas_size_t *N, blas_scomplex_t *A,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_scomplex_t *B,
       __LILA_BLAS_LAPACK_CONST blas_size_t *ldb, blas_float_t *W,
       blas_scomplex_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_float_t *rwork, blas_size_t *info);

extern "C" lapack_ret_t
zhegv_(__LILA_BLAS_LAPACK_CONST blas_size_t *itype,
       __LILA_BLAS_LAPACK_CONST char *jobz, __LILA_BLAS_LAPACK_CONST char *uplo,
       __LILA_BLAS_LAPACK_CONST blas_size_t *N, blas_complex_t *A,
       __LILA_BLAS_LAPACK_CONST blas_size_t *lda, blas_complex_t *B,
       __LILA_BLAS_LAPACK_CONST blas_size_t *ldb, blas_double_t *W,
       blas_complex_t *work, __LILA_BLAS_LAPACK_CONST blas_size_t *lwork,
       blas_double_t *rwork, blas_size_t *info);

#endif
