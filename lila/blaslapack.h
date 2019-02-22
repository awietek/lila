// Copyright 2018 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LILA_BLASLAPACK_H_
#define LILA_BLASLAPACK_H_

#include <cstdlib>
#include <complex>

#ifdef LILA_USE_MKL // Using the Intel MKL
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "blaslapack_types.h"
#include "mkl_types.h"
#include "mkl.h"
#define __LAPACK_ROUTINE_NAME(x) x

#else  // Using normal LAPACK  
#include "blaslapack_types.h"
#include "blaslapack_extern.h"
#define __LAPACK_ROUTINE_NAME(x) x##_
#endif

namespace lila { 
  namespace blaslapack {

    // using blas_size_t = int;
    // using blas_complex_t = lila::blaslapack::blas_complex_t;

    // Copy
    inline void copy(const blas_size_t* N, const blas_float_t* x, 
		     const blas_size_t* incx, blas_float_t* y, 
		     const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(scopy)(N, x, incx, y, incy); }
    inline void copy(const blas_size_t* N, const blas_double_t* x, 
		     const blas_size_t* incx, blas_double_t* y, 
		     const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(dcopy)(N, x, incx, y, incy); }
    inline void copy(const blas_size_t* N, const blas_scomplex_t* x, 
		     const blas_size_t* incx, blas_scomplex_t* y, 
		     const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(ccopy)(N, x, incx, y, incy); }
    inline void copy(const blas_size_t* N, const blas_complex_t* x, 
		     const blas_size_t* incx, blas_complex_t* y, 
		     const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(zcopy)(N, x, incx, y, incy); }

    // Axpy
    inline void axpy(const blas_size_t* N,const blas_float_t* alpha, 
		     const blas_float_t* x, const blas_size_t* incx, 
		     blas_float_t* y,const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(saxpy)(N, alpha, x, incx, y, incy); }
    inline void axpy(const blas_size_t* N,const blas_double_t* alpha, 
		     const blas_double_t* x, const blas_size_t* incx, 
		     blas_double_t* y,const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(daxpy)(N, alpha, x, incx, y, incy); }
    inline void axpy(const blas_size_t* N,const blas_scomplex_t* alpha, 
		     const blas_scomplex_t* x, const blas_size_t* incx, 
		     blas_scomplex_t* y,const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(caxpy)(N, alpha, x, incx, y, incy); }
    inline void axpy(const blas_size_t* N,const blas_complex_t* alpha, 
		     const blas_complex_t* x, const blas_size_t* incx, 
		     blas_complex_t* y,const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(zaxpy)(N, alpha, x, incx, y, incy); }

    // Scal
    inline void scal(const blas_size_t* N,const blas_float_t* alpha, 
		     blas_float_t* x, const blas_size_t* incx)
    { __LAPACK_ROUTINE_NAME(sscal)(N, alpha, x, incx); }
    inline void scal(const blas_size_t* N,const blas_double_t* alpha, 
		     blas_double_t* x, const blas_size_t* incx)
    { __LAPACK_ROUTINE_NAME(dscal)(N, alpha, x, incx); }
    inline void scal(const blas_size_t* N,const blas_scomplex_t* alpha, 
		     blas_scomplex_t* x, const blas_size_t* incx)
    { __LAPACK_ROUTINE_NAME(cscal)(N, alpha, x, incx); }
    inline void scal(const blas_size_t* N,const blas_complex_t* alpha, 
		     blas_complex_t* x, const blas_size_t* incx)
    { __LAPACK_ROUTINE_NAME(zscal)(N, alpha, x, incx); }


    
    // Dot
    inline float dot(const blas_size_t* N ,const blas_float_t* x,
		     const blas_size_t* incx, const blas_float_t* y,
		     const blas_size_t* incy)
    { return __LAPACK_ROUTINE_NAME(sdot)(N, x, incx, y, incy); }
    inline double dot(const blas_size_t* N ,const blas_double_t* x,
		      const blas_size_t* incx, const blas_double_t* y,
		      const blas_size_t* incy)
    { return __LAPACK_ROUTINE_NAME(ddot)(N, x, incx, y, incy); }

    inline blas_scomplex_t dot(const blas_size_t* N, const blas_scomplex_t* x,
    			       const blas_size_t* incx, const blas_scomplex_t* y,
    			       const blas_size_t* incy)
    { 
#ifdef LILA_USE_MKL
      blas_scomplex_t res;
      __LAPACK_ROUTINE_NAME(cdotc)(&res, N, x, incx, y, incy);
      return res;
#else
      return __LAPACK_ROUTINE_NAME(cdotc)(N, x, incx, y, incy); 
#endif
    }
    
    inline blas_complex_t dot(const blas_size_t* N, const blas_complex_t* x,
			      const blas_size_t* incx, const blas_complex_t* y,
			      const blas_size_t* incy)
    { 
#ifdef LILA_USE_MKL
      blas_complex_t res;
      __LAPACK_ROUTINE_NAME(zdotc)(&res, N, x, incx, y, incy);
      return res;
#else
      return __LAPACK_ROUTINE_NAME(zdotc)(N, x, incx, y, incy); 
#endif
    }
        
    // Gemv
    inline void gemv(const char* trans, const blas_size_t* m, 
		     const blas_size_t* n, const blas_float_t* alpha, 
		     const blas_float_t* A, const blas_size_t* dima,
		     const blas_float_t* x, const blas_size_t* incx, 
		     const blas_float_t* beta, blas_float_t* y, 
		     const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(sgemv)(trans, m, n, alpha, A, dima, x, incx, 
				   beta, y, incy); }
    inline void gemv(const char* trans, const blas_size_t* m, 
		     const blas_size_t* n, const blas_double_t* alpha, 
		     const blas_double_t* A, const blas_size_t* dima,
		     const blas_double_t* x, const blas_size_t* incx, 
		     const blas_double_t* beta, blas_double_t* y, 
		     const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(dgemv)(trans, m, n, alpha, A, dima, x, incx, 
				   beta, y, incy); }
    inline void gemv(const char* trans, const blas_size_t* m, 
		     const blas_size_t* n, const blas_scomplex_t* alpha, 
		     const blas_scomplex_t* A, const blas_size_t* dima, 
		     const  blas_scomplex_t* x, const blas_size_t* incx, 
		     const blas_scomplex_t* beta, blas_scomplex_t* y, 
		     const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(cgemv)(trans, m, n, alpha, A, dima, x, incx, 
				   beta, y, incy); }
    inline void gemv(const char* trans, const blas_size_t* m,
		     const blas_size_t* n, const blas_complex_t* alpha, 
		     const blas_complex_t* A, const blas_size_t* dima,
		     const blas_complex_t* x, const blas_size_t* incx, 
		     const blas_complex_t* beta, blas_complex_t* y, 
		     const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(zgemv)(trans, m, n, alpha, A, dima, x, incx, 
				   beta, y, incy); }
    
    // Gemm
    inline void gemm(const char* transa, const char* transb, 
		     const blas_size_t* m, const blas_size_t* n, 
		     const blas_size_t* k, const blas_float_t* alpha,
		     const blas_float_t* A, const blas_size_t* dima, 
		     const blas_float_t* B, const blas_size_t* dimb, 
		     const blas_float_t* beta, blas_float_t* C,
		     const blas_size_t* dimc)
    { __LAPACK_ROUTINE_NAME(sgemm)(transa, transb, m, n, k, alpha, A, dima, 
				   B, dimb, beta, C, dimc); }
    inline void gemm(const char* transa, const char* transb, 
		     const blas_size_t* m, const blas_size_t* n, 
		     const blas_size_t* k, const blas_double_t* alpha,
		     const blas_double_t* A, const blas_size_t* dima, 
		     const blas_double_t* B, const blas_size_t* dimb, 
		     const blas_double_t* beta, blas_double_t* C,
		     const blas_size_t* dimc)
    { __LAPACK_ROUTINE_NAME(dgemm)(transa, transb, m, n, k, alpha, A, dima, 
				   B, dimb, beta, C, dimc); }
    inline void gemm(const char* transa, const char* transb, 
		     const blas_size_t* m, const blas_size_t* n, 
		     const blas_size_t* k, const blas_scomplex_t* alpha,
		     const blas_scomplex_t* A, const blas_size_t* dima, 
		     const blas_scomplex_t* B, const blas_size_t* dimb, 
		     const blas_scomplex_t* beta, blas_scomplex_t* C,
		     const blas_size_t* dimc)  
    { __LAPACK_ROUTINE_NAME(cgemm)(transa, transb, m, n, k, alpha, A, dima, 
				   B, dimb, beta, C, dimc); }
    inline void gemm(const char* transa, const char* transb, 
		     const blas_size_t* m, const blas_size_t* n, 
		     const blas_size_t* k, const blas_complex_t* alpha,
		     const blas_complex_t* A, const blas_size_t* dima, 
		     const blas_complex_t* B, const blas_size_t* dimb, 
		     const blas_complex_t* beta, blas_complex_t* C,
		     const blas_size_t* dimc)
    { __LAPACK_ROUTINE_NAME(zgemm)(transa, transb, m, n, k, alpha, A, dima, 
				   B, dimb, beta, C, dimc); }

    //////////////////////////
    // Linear Solve
    // Gesv
    inline void gesv(const blas_size_t* n, const blas_size_t* n_rhs, 
		     blas_float_t* A, const blas_size_t* lda, 
		     blas_size_t* ipiv, blas_float_t* B,
		     const blas_size_t* ldb, int* info)
    { __LAPACK_ROUTINE_NAME(sgesv)(n, n_rhs, A, lda, ipiv, B, ldb, info); }
    inline void gesv(const blas_size_t* n, const blas_size_t* n_rhs, 
		     blas_double_t* A, const blas_size_t* lda, 
		     blas_size_t* ipiv, blas_double_t* B,
		     const blas_size_t* ldb, int* info)
    { __LAPACK_ROUTINE_NAME(dgesv)(n, n_rhs, A, lda, ipiv, B, ldb, info); }
    inline void gesv(const blas_size_t* n, const blas_size_t* n_rhs, 
		     blas_scomplex_t* A, const blas_size_t* lda, 
		     blas_size_t* ipiv, blas_scomplex_t* B,
		     const blas_size_t* ldb, int* info)
    { __LAPACK_ROUTINE_NAME(cgesv)(n, n_rhs, A, lda, ipiv, B, ldb, info); }
    inline void gesv(const blas_size_t* n, const blas_size_t* n_rhs, 
		     blas_complex_t* A, const blas_size_t* lda, 
		     blas_size_t* ipiv, blas_complex_t* B,
		     const blas_size_t* ldb, int* info)
    { __LAPACK_ROUTINE_NAME(zgesv)(n, n_rhs, A, lda, ipiv, B, ldb, info); }

    //////////////////////////
    // LU Decomposition
    // Getrf (performs LU Decomposition)
    inline void getrf(const blas_size_t* M, const blas_size_t* N, 
		      blas_float_t* A, const blas_size_t* lda,
		      blas_size_t* ipiv, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(sgetrf)(M, N, A, lda, ipiv, info); }
    inline void getrf(const blas_size_t* M, const blas_size_t* N, 
		      blas_double_t* A, const blas_size_t* lda,
		      blas_size_t* ipiv, blas_size_t* info)    
    { __LAPACK_ROUTINE_NAME(dgetrf)(M, N, A, lda, ipiv, info); }
    inline void getrf(const blas_size_t* M, const blas_size_t* N, 
		      blas_scomplex_t* A, const blas_size_t* lda,
		      blas_size_t* ipiv, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(cgetrf)(M, N, A, lda, ipiv, info); }
    inline void getrf(const blas_size_t* M, const blas_size_t* N, 
		      blas_complex_t* A, const blas_size_t* lda,
		      blas_size_t* ipiv, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(zgetrf)(M, N, A, lda, ipiv, info); }

    // Getrs (solves system of equations from LU Decomposition)
    inline void getrs(const char* trans, const blas_size_t* n, 
		      const blas_size_t* n_rhs, const blas_float_t* A,
		      const blas_size_t* lda, const blas_size_t* ipiv,
		      blas_float_t* B, const blas_size_t* ldb,
		      blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(sgetrs)(trans, n, n_rhs, A, lda, ipiv, B, 
				    ldb, info); }
    inline void getrs(const char* trans, const blas_size_t* n, 
		      const blas_size_t* n_rhs, const blas_double_t* A,
		      const blas_size_t* lda, const blas_size_t* ipiv,
		      blas_double_t* B, const blas_size_t* ldb,
		      blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(dgetrs)(trans, n, n_rhs, A, lda, ipiv, B, 
				    ldb, info); }
    inline void getrs(const char* trans, const blas_size_t* n, 
		      const blas_size_t* n_rhs, const blas_scomplex_t* A,
		      const blas_size_t* lda, const blas_size_t* ipiv,
		      blas_scomplex_t* B, const blas_size_t* ldb,
		      blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(cgetrs)(trans, n, n_rhs, A, lda, ipiv, B, 
				    ldb, info); }
    inline void getrs(const char* trans, const blas_size_t* n, 
		      const blas_size_t* n_rhs, const blas_complex_t* A,
		      const blas_size_t* lda, const blas_size_t* ipiv,
		      blas_complex_t* B, const blas_size_t* ldb,
		      blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(zgetrs)(trans, n, n_rhs, A, lda, ipiv, B, 
				    ldb, info); }

    // Getri (computes inverse from LU Decomposition)
    inline void getri(const blas_size_t* N, blas_float_t* A, 
		      const blas_size_t* lda, const blas_size_t* IPIV, 
		      float* WORK, const blas_size_t* lwork, 
		      blas_size_t* INFO)    
    { __LAPACK_ROUTINE_NAME(sgetri)(N, A, lda, IPIV, WORK, lwork, INFO); }
    inline void getri(const blas_size_t* N, blas_double_t* A, 
		      const blas_size_t* lda, const blas_size_t* IPIV, 
		      blas_double_t* WORK, const blas_size_t* lwork, 
		      blas_size_t* INFO)
    { __LAPACK_ROUTINE_NAME(dgetri)(N, A, lda, IPIV, WORK, lwork, INFO); }
    inline void getri(const blas_size_t* N, blas_scomplex_t* A, 
		      const blas_size_t* lda, const blas_size_t* IPIV, 
		      blas_scomplex_t* WORK, const blas_size_t* lwork, 
		      blas_size_t* INFO)
    { __LAPACK_ROUTINE_NAME(cgetri)(N, A, lda, IPIV, WORK, lwork, INFO); }
    inline void getri(const blas_size_t* N, blas_complex_t* A, 
		      const blas_size_t* lda, const blas_size_t* IPIV, 
		      blas_complex_t* WORK, const blas_size_t* lwork, 
		      blas_size_t* INFO)
    { __LAPACK_ROUTINE_NAME(zgetri)(N, A, lda, IPIV, WORK, lwork, INFO); }


    //////////////////////////
    // QR Decomposition
    // Geqrf  (does the QR Decomposition)
    inline void geqrf(const blas_size_t* m, const blas_size_t* n,	
		      blas_float_t* A, const blas_size_t* lda,
		      blas_float_t* tau, blas_float_t* work, 
		      const blas_size_t* lwork, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(sgeqrf)(m, n, A, lda, tau, work, lwork, info); }
    inline void geqrf(const blas_size_t* m, const blas_size_t* n,	
		      blas_double_t* A, const blas_size_t* lda,
		      blas_double_t* tau, blas_double_t* work, 
		      const blas_size_t* lwork, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(dgeqrf)(m, n, A, lda, tau, work, lwork, info); }
    inline void geqrf(const blas_size_t* m, const blas_size_t* n,	
		      blas_scomplex_t* A, const blas_size_t* lda,
		      blas_scomplex_t* tau, blas_scomplex_t* work, 
		      const blas_size_t* lwork, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(cgeqrf)(m, n, A, lda, tau, work, lwork, info); }
    inline void geqrf(const blas_size_t* m, const blas_size_t* n,	
		      blas_complex_t* A, const blas_size_t* lda,
		      blas_complex_t* tau, blas_complex_t* work, 
		      const blas_size_t* lwork, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(zgeqrf)(m, n, A, lda, tau, work, lwork, info); }

    // Orgqr / Ungqr (retrieves the matrix Q of the QR Decomposition)
    // (Note: use orgqr also for the complex routines xungqr)
    inline void orgqr(const blas_size_t* m, const blas_size_t* n,
		      const blas_size_t* k, blas_float_t* A, 
		      const blas_size_t* lda, const blas_float_t* tau, 
		      blas_float_t* work, const blas_size_t* lwork, 
		      blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(sorgqr)(m, n, k, A, lda, tau, work, lwork, info); }
    inline void orgqr(const blas_size_t* m, const blas_size_t* n,
		      const blas_size_t* k, blas_double_t* A, 
		      const blas_size_t* lda, const blas_double_t* tau, 
		      blas_double_t* work, const blas_size_t* lwork, 
		      blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(dorgqr)(m, n, k, A, lda, tau, work, lwork, info); }
    inline void orgqr(const blas_size_t* m, const blas_size_t* n,
		      const blas_size_t* k, blas_scomplex_t* A, 
		      const blas_size_t* lda, const blas_scomplex_t* tau, 
		      blas_scomplex_t* work, const blas_size_t* lwork, 
		      blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(cungqr)(m, n, k, A, lda, tau, work, lwork, info); }
    inline void orgqr(const blas_size_t* m, const blas_size_t* n,
		      const blas_size_t* k, blas_complex_t* A, 
		      const blas_size_t* lda, const blas_complex_t* tau, 
		      blas_complex_t* work, const blas_size_t* lwork, 
		      blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(zungqr)(m, n, k, A, lda, tau, work, lwork, info); }
      
    //////////////////////////
    // Eigenvalues

    // symmetric/hermitian
    inline void syev(const char* jobz, const char* uplo, 
		     const blas_size_t* n, blas_float_t* a, 
		     const blas_size_t* lda, blas_float_t* w, 
		     blas_float_t* work, const blas_size_t* lwork, 
		     blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(ssyev)(jobz, uplo, n, a, lda, 
				   w, work, lwork, info); }
    
    inline void syev(const char* jobz, const char* uplo, 
		     const blas_size_t* n, blas_double_t* a, 
		     const blas_size_t* lda, blas_double_t* w, 
		     blas_double_t* work, const blas_size_t* lwork, 
		     blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(dsyev)(jobz, uplo, n, a, lda, 
				   w, work, lwork, info); }

    inline void heev(const char* jobz, const char* uplo, 
		     const blas_size_t* n, blas_scomplex_t* a, 
		     const blas_size_t* lda, blas_float_t* w, 
		     blas_scomplex_t* work, const blas_size_t* lwork, 
		     blas_float_t* rwork, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(cheev)(jobz, uplo, n, a, lda, 
				   w, work, lwork, rwork, info); }
      
    inline void heev(const char* jobz, const char* uplo, 
		     const blas_size_t* n, blas_complex_t* a, 
		     const blas_size_t* lda, blas_double_t* w, 
		     blas_complex_t* work, const blas_size_t* lwork, 
		     blas_double_t* rwork, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(zheev)(jobz, uplo, n, a, lda, 
				   w, work, lwork, rwork, info); }

    // (not actual lapack routines, create extra rwork)
    inline void syev(const char* jobz, const char* uplo, 
		     const blas_size_t* n, blas_scomplex_t* a, 
		     const blas_size_t* lda, blas_float_t* w, 
		     blas_scomplex_t* work, const blas_size_t* lwork,
		     blas_size_t* info)
    { 
      std::vector<blas_float_t> rwork(3*(*n)-2); 
      heev(jobz, uplo, n, a, lda, w, work, lwork, rwork.data(), info);
    }
      
    inline void syev(const char* jobz, const char* uplo, 
		     const blas_size_t* n, blas_complex_t* a, 
		     const blas_size_t* lda, blas_double_t* w, 
		     blas_complex_t* work, const blas_size_t* lwork,
		     blas_size_t* info)
    { 
      std::vector<blas_double_t> rwork(3*(*n)-2); 
      heev(jobz, uplo, n, a, lda, w, work, lwork, rwork.data(), info);
    }

    // generic
    inline void geev(const char* jobvl, const char* jobvr,
		     const blas_size_t* n, blas_float_t* a, 
		     const blas_size_t* lda, blas_float_t* wr, 
		     blas_float_t* wi, blas_float_t* vl,
		     const blas_size_t* ldvl, blas_float_t* vr,
		     const blas_size_t* ldvr, blas_float_t* work,
		     const blas_size_t* lwork, blas_size_t* info)
    {  __LAPACK_ROUTINE_NAME(sgeev)(jobvl, jobvr, n, a, lda, wr, wi, vl, 
				   ldvl, vr, ldvr, work, lwork, info); }
    inline void geev(const char* jobvl, const char* jobvr,
		     const blas_size_t* n, blas_double_t* a, 
		     const blas_size_t* lda, blas_double_t* wr, 
		     blas_double_t* wi, blas_double_t* vl,
		     const blas_size_t* ldvl, blas_double_t* vr,
		     const blas_size_t* ldvr, blas_double_t* work,
		     const blas_size_t* lwork, blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(dgeev)(jobvl, jobvr, n, a, lda, wr, wi, vl, 
				   ldvl, vr, ldvr, work, lwork, info); }
    
    // (not actual lapack routines, convert real/imag vectors to complex)
    inline void geev(const char* jobvl, const char* jobvr,
		     const blas_size_t* n, blas_float_t* a, 
		     const blas_size_t* lda, blas_scomplex_t* w,  
		     blas_float_t* vl, const blas_size_t* ldvl, 
		     blas_float_t* vr, const blas_size_t* ldvr, 
		     blas_float_t* work, const blas_size_t* lwork, 
		     blas_size_t* info)
    { 
      std::vector<blas_float_t> wr(*n);
      std::vector<blas_float_t> wi(*n);
      geev(jobvl, jobvr, n, a, lda, wr.data(), wi.data(), vl, ldvl, vr, 
	   ldvr, work, lwork, info); 
      for (blas_size_t i=0; i<*n; ++i)
	w[i] = blas_scomplex_t({wr[i], wi[i]});      
    }

    inline void geev(const char* jobvl, const char* jobvr,
		     const blas_size_t* n, blas_double_t* a, 
		     const blas_size_t* lda, blas_complex_t* w, 
		     blas_double_t* vl, const blas_size_t* ldvl, 
		     blas_double_t* vr, const blas_size_t* ldvr, 
		     blas_double_t* work, const blas_size_t* lwork, 
		     blas_size_t* info)
    {       
      std::vector<blas_double_t> wr(*n);
      std::vector<blas_double_t> wi(*n);
      geev(jobvl, jobvr, n, a, lda, wr.data(), wi.data(), vl, ldvl, vr, 
	   ldvr, work, lwork, info); 
      for (blas_size_t i=0; i<*n; ++i)
	w[i] = blas_complex_t({wr[i], wi[i]});  
    }

    inline void geev(const char* jobvl, const char* jobvr,
		     const blas_size_t* n, blas_scomplex_t* a, 
		     const blas_size_t* lda, blas_scomplex_t* w, 
		     blas_scomplex_t* vl, const blas_size_t* ldvl, 
		     blas_scomplex_t* vr, const blas_size_t* ldvr, 
		     blas_scomplex_t* work, const blas_size_t* lwork, 
		     blas_size_t* info)
    { 
      std::vector<blas_float_t> rwork(*n * 2);
      __LAPACK_ROUTINE_NAME(cgeev)(jobvl, jobvr, n, a, lda, w, vl, ldvl, 
				   vr, ldvr, work, lwork, rwork.data(),info); 
    }
    
    inline void geev(const char* jobvl, const char* jobvr,
		     const blas_size_t* n, blas_complex_t* a, 
		     const blas_size_t* lda, blas_complex_t* w, 
		     blas_complex_t* vl, const blas_size_t* ldvl, 
		     blas_complex_t* vr, const blas_size_t* ldvr, 
		     blas_complex_t* work, const blas_size_t* lwork, 
		     blas_size_t* info)
    {
      std::vector<blas_double_t> rwork(*n * 2);  
      __LAPACK_ROUTINE_NAME(zgeev)(jobvl, jobvr, n, a, lda, w, vl, ldvl,
				   vr, ldvr, work, lwork, rwork.data(), info); 
    }


    // Real symmetric tridiagonal eigensolvers
    inline void stev(const char* jobz, const blas_size_t* N, 
		     blas_float_t* D, blas_float_t* E, blas_float_t* Z,
		     const blas_size_t* ldz, blas_float_t* work, 
		     blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(sstev)(jobz, N, D, E, Z, ldz, work, info); }

    inline void stev(const char* jobz, const blas_size_t* N, 
		     blas_double_t* D, blas_double_t* E, blas_double_t* Z,
		     const blas_size_t* ldz, blas_double_t* work, 
		     blas_size_t* info)
    { __LAPACK_ROUTINE_NAME(dstev)(jobz, N, D, E, Z, ldz, work, info); }
    
  }  // namespace lila
}  // namespace blaslapack

#endif
