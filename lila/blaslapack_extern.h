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

#ifndef LILA_BLASLAPACK_EXTERN_H_
#define LILA_BLASLAPACK_EXTERN_H_

#include "blaslapack_types.h"

#ifndef LILA_USE_MKL

namespace lila { 
  namespace blaslapack {
    
    // Copy
    extern "C" void scopy_(const blas_size_t* N, const blas_float_t* x, 
			   const blas_size_t* incx, blas_float_t* y, 
			   const blas_size_t* incy);
    extern "C" void dcopy_(const blas_size_t* N, const blas_double_t* x, 
			   const blas_size_t* incx, blas_double_t* y, 
			   const blas_size_t* incy);
    extern "C" void ccopy_(const blas_size_t* N, const blas_scomplex_t* x, 
			   const blas_size_t* incx, blas_scomplex_t* y, 
			   const blas_size_t* incy);
    extern "C" void zcopy_(const blas_size_t* N, const blas_complex_t* x, 
			   const blas_size_t* incx, blas_complex_t* y, 
			   const blas_size_t* incy);

    // Axpy
    extern "C" void saxpy_(const blas_size_t* N,const blas_float_t* alpha, 
			   const blas_float_t* x, const blas_size_t* incx, 
			   blas_float_t* y,const blas_size_t* incy);
    extern "C" void daxpy_(const blas_size_t* N,const blas_double_t* alpha, 
			   const blas_double_t* x, const blas_size_t* incx, 
			   blas_double_t* y,const blas_size_t* incy);
    extern "C" void caxpy_(const blas_size_t* N,const blas_scomplex_t* alpha, 
			   const blas_scomplex_t* x, const blas_size_t* incx, 
			   blas_scomplex_t* y,const blas_size_t* incy);
    extern "C" void zaxpy_(const blas_size_t* N,const blas_complex_t* alpha, 
			   const blas_complex_t* x, const blas_size_t* incx, 
			   blas_complex_t* y,const blas_size_t* incy);

    // Scal
    extern "C" void sscal_(const blas_size_t* N,const blas_float_t* alpha, 
			   blas_float_t* x, const blas_size_t* incx);
    extern "C" void dscal_(const blas_size_t* N,const blas_double_t* alpha, 
			   blas_double_t* x, const blas_size_t* incx);
    extern "C" void cscal_(const blas_size_t* N,const blas_scomplex_t* alpha, 
			   blas_scomplex_t* x, const blas_size_t* incx);
    extern "C" void zscal_(const blas_size_t* N,const blas_complex_t* alpha, 
			   blas_complex_t* x, const blas_size_t* incx);

    // Dot
    extern "C" float sdot_(const blas_size_t* N ,const blas_float_t* x,
			   const blas_size_t* incx, const blas_float_t* y,
			   const blas_size_t* incy);
    extern "C" double ddot_(const blas_size_t* N ,const blas_double_t* x,
			   const blas_size_t* incx, const blas_double_t* y,
			   const blas_size_t* incy);
    extern "C" blas_scomplex_t cdotc_(const blas_size_t* N,
				      const blas_scomplex_t* x,
				      const blas_size_t* incx, 
				      const blas_scomplex_t* y,
				      const blas_size_t* incy);
    extern "C" blas_complex_t zdotc_(const blas_size_t* N,
				     const blas_complex_t* x,
				     const blas_size_t* incx, 
				     const blas_complex_t* y,
				     const blas_size_t* incy);

    // Gemv
    extern "C" void sgemv_(const char* trans, const blas_size_t* m, 
			   const blas_size_t* n, const blas_float_t* alpha, 
			   const blas_float_t* A, const blas_size_t* dima,
			   const blas_float_t* x, const blas_size_t* incx, 
			   const blas_float_t* beta, blas_float_t* y, 
			   const blas_size_t* incy);
    extern "C" void dgemv_(const char* trans, const blas_size_t* m, 
			   const blas_size_t* n, const blas_double_t* alpha, 
			   const blas_double_t* A, const blas_size_t* dima,
			   const blas_double_t* x, const blas_size_t* incx, 
			   const blas_double_t* beta, blas_double_t* y, 
			   const blas_size_t* incy);
    extern "C" void cgemv_(const char* trans, const blas_size_t* m, 
			   const blas_size_t* n, const blas_scomplex_t* alpha, 
			   const blas_scomplex_t* A, const blas_size_t* dima, 
			   const  blas_scomplex_t* x, const blas_size_t* incx, 
			   const blas_scomplex_t* beta, blas_scomplex_t* y, 
			   const blas_size_t* incy);
    extern "C" void zgemv_(const char* trans, const blas_size_t* m,
			   const blas_size_t* n, const blas_complex_t* alpha, 
			   const blas_complex_t* A, const blas_size_t* dima,
			   const blas_complex_t* x, const blas_size_t* incx, 
			   const blas_complex_t* beta, blas_complex_t* y, 
			   const blas_size_t* incy);
        
    // Gemm
    extern "C" void sgemm_(const char* transa, const char* transb, 
			   const blas_size_t* m, const blas_size_t* k, 
			   const blas_size_t* n, const blas_float_t* alpha,
			   const blas_float_t* A, const blas_size_t* dima, 
			   const blas_float_t* B, const blas_size_t* dimb, 
			   const blas_float_t* beta, blas_float_t* C,
			   const blas_size_t* dimc);
    extern "C" void dgemm_(const char* transa, const char* transb, 
			   const blas_size_t* m, const blas_size_t* k, 
			   const blas_size_t* n, const blas_double_t* alpha,
			   const blas_double_t* A, const blas_size_t* dima, 
			   const blas_double_t* B, const blas_size_t* dimb, 
			   const blas_double_t* beta, blas_double_t* C,
			   const blas_size_t* dimc);
    extern "C" void cgemm_(const char* transa, const char* transb, 
			   const blas_size_t* m, const blas_size_t* k, 
			   const blas_size_t* n, const blas_scomplex_t* alpha,
			   const blas_scomplex_t* A, const blas_size_t* dima, 
			   const blas_scomplex_t* B, const blas_size_t* dimb, 
			   const blas_scomplex_t* beta, blas_scomplex_t* C,
			   const blas_size_t* dimc);  
    extern "C" void zgemm_(const char* transa, const char* transb, 
			   const blas_size_t* m, const blas_size_t* k, 
			   const blas_size_t* n, const blas_complex_t* alpha,
			   const blas_complex_t* A, const blas_size_t* dima, 
			   const blas_complex_t* B, const blas_size_t* dimb, 
			   const blas_complex_t* beta, blas_complex_t* C,
			   const blas_size_t* dimc);

    //////////////////////////
    // Linear Solve
    // Gesv
    extern "C" void sgesv_(const blas_size_t* n, const blas_size_t* n_rhs, 
			   blas_float_t* A, const blas_size_t* lda, 
			   blas_size_t* ipiv, blas_float_t* B,
			   const blas_size_t* ldb, int* info);
    extern "C" void dgesv_(const blas_size_t* n, const blas_size_t* n_rhs, 
			   blas_double_t* A, const blas_size_t* lda, 
			   blas_size_t* ipiv, blas_double_t* B,
			   const blas_size_t* ldb, int* info);
    extern "C" void cgesv_(const blas_size_t* n, const blas_size_t* n_rhs, 
			   blas_scomplex_t* A, const blas_size_t* lda, 
			   blas_size_t* ipiv, blas_scomplex_t* B,
			   const blas_size_t* ldb, int* info);
    extern "C" void zgesv_(const blas_size_t* n, const blas_size_t* n_rhs, 
			   blas_complex_t* A, const blas_size_t* lda, 
			   blas_size_t* ipiv, blas_complex_t* B,
			   const blas_size_t* ldb, int* info);

    //////////////////////////
    // LU Decomposition
    // Getrf (performs LU Decomposition)
    extern "C" void sgetrf_(const blas_size_t* M, const blas_size_t* N, 
			    blas_float_t* A, const blas_size_t* lda,
			    blas_size_t* ipiv, blas_size_t* info);
    extern "C" void dgetrf_(const blas_size_t* M, const blas_size_t* N, 
			    blas_double_t* A, const blas_size_t* lda,
			    blas_size_t* ipiv, blas_size_t* info);
    extern "C" void cgetrf_(const blas_size_t* M, const blas_size_t* N, 
			    blas_scomplex_t* A, const blas_size_t* lda,
			    blas_size_t* ipiv, blas_size_t* info);
    extern "C" void zgetrf_(const blas_size_t* M, const blas_size_t* N, 
			    blas_complex_t* A, const blas_size_t* lda,
			    blas_size_t* ipiv, blas_size_t* info);
    
    // Getrs (solves system of equations from LU Decomposition)
    extern "C" void sgetrs_(const char* trans, const blas_size_t* n, 
			    const blas_size_t* n_rhs, const blas_float_t* A,
			    const blas_size_t* lda, const blas_size_t* ipiv,
			    blas_float_t* B, const blas_size_t* ldb,
			    blas_size_t* info);
    extern "C" void dgetrs_(const char* trans, const blas_size_t* n, 
			    const blas_size_t* n_rhs, const blas_double_t* A,
			    const blas_size_t* lda, const blas_size_t* ipiv,
			    blas_double_t* B, const blas_size_t* ldb,
			    blas_size_t* info);
    extern "C" void cgetrs_(const char* trans, const blas_size_t* n, 
			    const blas_size_t* n_rhs, const blas_scomplex_t* A,
			    const blas_size_t* lda, const blas_size_t* ipiv,
			    blas_scomplex_t* B, const blas_size_t* ldb,
			    blas_size_t* info);
    extern "C" void zgetrs_(const char* trans, const blas_size_t* n, 
			    const blas_size_t* n_rhs, const blas_complex_t* A,
			    const blas_size_t* lda, const blas_size_t* ipiv,
			    blas_complex_t* B, const blas_size_t* ldb,
			    blas_size_t* info);

    // Getri (computes inverse from LU Decomposition)
    extern "C" void sgetri_(const blas_size_t* N, blas_float_t* A, 
			    const blas_size_t* lda, const blas_size_t* IPIV, 
			    blas_float_t* WORK, const blas_size_t* lwork, 
			    blas_size_t* INFO);
    extern "C" void dgetri_(const blas_size_t* N, blas_double_t* A, 
			    const blas_size_t* lda, const blas_size_t* IPIV, 
			    blas_double_t* WORK, const blas_size_t* lwork, 
			    blas_size_t* INFO);
    extern "C" void cgetri_(const blas_size_t* N, blas_scomplex_t* A, 
			    const blas_size_t* lda, const blas_size_t* IPIV, 
			    blas_scomplex_t* WORK, const blas_size_t* lwork, 
			    blas_size_t* INFO);
    extern "C" void zgetri_(const blas_size_t* N, blas_complex_t* A, 
			    const blas_size_t* lda, const blas_size_t* IPIV, 
			    blas_complex_t* WORK, const blas_size_t* lwork, 
			    blas_size_t* INFO);

    //////////////////////////
    // QR Decomposition
    // Geqrf (does the QR Decomposition)
    extern "C" void sgeqrf_(const blas_size_t* m, const blas_size_t* n,	
			    blas_float_t* A, const blas_size_t* lda,
			    blas_float_t* tau, blas_float_t* work, 
			    const blas_size_t* lwork, blas_size_t* info);
    extern "C" void dgeqrf_(const blas_size_t* m, const blas_size_t* n,	
			    blas_double_t* A, const blas_size_t* lda,
			    blas_double_t* tau, blas_double_t* work, 
			    const blas_size_t* lwork, blas_size_t* info);
    extern "C" void cgeqrf_(const blas_size_t* m, const blas_size_t* n,	
			    blas_scomplex_t* A, const blas_size_t* lda,
			    blas_scomplex_t* tau, blas_scomplex_t* work, 
			    const blas_size_t* lwork, blas_size_t* info);
    extern "C" void zgeqrf_(const blas_size_t* m, const blas_size_t* n,	
			    blas_complex_t* A, const blas_size_t* lda,
			    blas_complex_t* tau, blas_complex_t* work, 
			    const blas_size_t* lwork, blas_size_t* info);

    // Orgqr / Ungqr (retrieves the matrix Q of the QR Decomposition)
    extern "C" void sorgqr_(const blas_size_t* m, const blas_size_t* n,
			    const blas_size_t* k, blas_float_t* A, 
			    const blas_size_t* lda, const blas_float_t* tau, 
			    blas_float_t* work, const blas_size_t* lwork, 
			    blas_size_t* info);
    extern "C" void dorgqr_(const blas_size_t* m, const blas_size_t* n,
			    const blas_size_t* k, blas_double_t* A, 
			    const blas_size_t* lda, const blas_double_t* tau, 
			    blas_double_t* work, const blas_size_t* lwork, 
			    blas_size_t* info);
    extern "C" void cungqr_(const blas_size_t* m, const blas_size_t* n,
			    const blas_size_t* k, blas_scomplex_t* A, 
			    const blas_size_t* lda, const blas_scomplex_t* tau, 
			    blas_scomplex_t* work, const blas_size_t* lwork, 
			    blas_size_t* info);
    extern "C" void zungqr_(const blas_size_t* m, const blas_size_t* n,
			    const blas_size_t* k, blas_complex_t* A, 
			    const blas_size_t* lda, const blas_complex_t* tau, 
			    blas_complex_t* work, const blas_size_t* lwork, 
			    blas_size_t* info);

    //////////////////////////
    // Eigenvalues

    // symmetric/hermitian
    extern "C" void ssyev_(const char* jobz, const char* uplo, 
			   const blas_size_t* n, blas_float_t* a, 
			   const blas_size_t* lda, blas_float_t* w, 
			   blas_float_t* work, const blas_size_t* lwork, 
			   blas_size_t* info);
    
    extern "C" void dsyev_(const char* jobz, const char* uplo, 
			   const blas_size_t* n, blas_double_t* a, 
			   const blas_size_t* lda, blas_double_t* w, 
			   blas_double_t* work, const blas_size_t* lwork, 
			   blas_size_t* info);

    extern "C" void cheev_(const char* jobz, const char* uplo, 
			   const blas_size_t* n, blas_scomplex_t* a, 
			   const blas_size_t* lda, blas_float_t* w, 
			   blas_scomplex_t* work, const blas_size_t* lwork, 
			   blas_float_t* rwork, blas_size_t* info);
    
    extern "C" void zheev_(const char* jobz, const char* uplo, 
			   const blas_size_t* n, blas_complex_t* a, 
			   const blas_size_t* lda, blas_double_t* w, 
			   blas_complex_t* work, const blas_size_t* lwork, 
			   blas_double_t* rwork, blas_size_t* info);

    // generic
    extern "C" void sgeev_(const char* jobvl, const char* jobvr,
			   const blas_size_t* n, blas_float_t* a, 
			   const blas_size_t* lda, blas_float_t* wr, 
			   blas_float_t* wi, blas_float_t* vl,
			   const blas_size_t* ldvl, blas_float_t* vr,
			   const blas_size_t* ldvr, blas_float_t* work,
			   const blas_size_t* lwork, blas_size_t* info);
    extern "C" void dgeev_(const char* jobvl, const char* jobvr,
			   const blas_size_t* n, blas_double_t* a, 
			   const blas_size_t* lda, blas_double_t* wr, 
			   blas_double_t* wi, blas_double_t* vl,
			   const blas_size_t* ldvl, blas_double_t* vr,
			   const blas_size_t* ldvr, blas_double_t* work,
			   const blas_size_t* lwork, blas_size_t* info);
    extern "C" void cgeev_(const char* jobvl, const char* jobvr,
			   const blas_size_t* n, blas_scomplex_t* a, 
			   const blas_size_t* lda, blas_scomplex_t* w, 
			   blas_scomplex_t* vl, const blas_size_t* ldvl, 
			   blas_scomplex_t* vr, const blas_size_t* ldvr, 
			   blas_scomplex_t* work, const blas_size_t* lwork, 
			   blas_float_t* rwork, blas_size_t* info);
    extern "C" void zgeev_(const char* jobvl, const char* jobvr,
			   const blas_size_t* n, blas_complex_t* a, 
			   const blas_size_t* lda, blas_complex_t* w, 
			   blas_complex_t* vl, const blas_size_t* ldvl, 
			   blas_complex_t* vr, const blas_size_t* ldvr, 
			   blas_complex_t* work, const blas_size_t* lwork, 
			   blas_double_t* rwork, blas_size_t* info);
    
    // Real symmetric tridiagonal eigensolvers
    extern "C" void sstev_(const char* jobz, const blas_size_t* N, 
			   blas_float_t* D, blas_float_t* E, blas_float_t* Z,
			   const blas_size_t* ldz, blas_float_t* work, 
			   blas_size_t* info);

    extern "C" void dstev_(const char* jobz, const blas_size_t* N, 
			   blas_double_t* D, blas_double_t* E, blas_double_t* Z,
			   const blas_size_t* ldz, blas_double_t* work, 
			   blas_size_t* info);

  }  // namespace blaslapack
}  // namespace lila

#endif

#endif
