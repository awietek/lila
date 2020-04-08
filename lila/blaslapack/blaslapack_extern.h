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

using blas_size_t = lila::blaslapack::blas_size_t;
using blas_float_t = lila::blaslapack::blas_float_t;
using blas_double_t = lila::blaslapack::blas_double_t;
using blas_scomplex_t = lila::blaslapack::blas_scomplex_t;
using blas_complex_t = lila::blaslapack::blas_complex_t;
using lapack_ret_t = lila::blaslapack::lapack_ret_t;

#if defined(LILA_USE_LAPACK) or defined(LILA_USE_ACCELERATE)



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
#ifdef LILA_USE_ACCELERATE
extern "C" void cdotc_(blas_scomplex_t* ret, const blas_size_t* N,
		       const blas_scomplex_t* x, const blas_size_t* incx, 
		       const blas_scomplex_t* y, const blas_size_t* incy);
extern "C" void zdotc_(blas_complex_t* ret, const blas_size_t* N,
		       const blas_complex_t* x, const blas_size_t* incx, 
		       const blas_complex_t* y, const blas_size_t* incy);
#else
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
#endif

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
#endif
    

#ifdef LILA_USE_LAPACK
//////////////////////////
// Linear Solve
// Gesv
extern "C" lapack_ret_t
sgesv_(LAPACK_CONST blas_size_t* n, LAPACK_CONST blas_size_t* n_rhs, 
       blas_float_t* A, LAPACK_CONST blas_size_t* lda, 
       blas_size_t* ipiv, blas_float_t* B,
       LAPACK_CONST blas_size_t* ldb, int* info);
extern "C" lapack_ret_t
dgesv_(LAPACK_CONST blas_size_t* n, LAPACK_CONST blas_size_t* n_rhs, 
       blas_double_t* A, LAPACK_CONST blas_size_t* lda, 
       blas_size_t* ipiv, blas_double_t* B,
       LAPACK_CONST blas_size_t* ldb, int* info);
extern "C" lapack_ret_t
cgesv_(LAPACK_CONST blas_size_t* n, LAPACK_CONST blas_size_t* n_rhs, 
       blas_scomplex_t* A, LAPACK_CONST blas_size_t* lda, 
       blas_size_t* ipiv, blas_scomplex_t* B,
       LAPACK_CONST blas_size_t* ldb, int* info);
extern "C" lapack_ret_t
zgesv_(LAPACK_CONST blas_size_t* n, LAPACK_CONST blas_size_t* n_rhs, 
       blas_complex_t* A, LAPACK_CONST blas_size_t* lda, 
       blas_size_t* ipiv, blas_complex_t* B,
       LAPACK_CONST blas_size_t* ldb, int* info);

//////////////////////////
// LU Decomposition
// Getrf (performs LU Decomposition)
extern "C" lapack_ret_t
sgetrf_(LAPACK_CONST blas_size_t* M, LAPACK_CONST blas_size_t* N, 
	blas_float_t* A, LAPACK_CONST blas_size_t* lda,
	blas_size_t* ipiv, blas_size_t* info);
extern "C" lapack_ret_t
dgetrf_(LAPACK_CONST blas_size_t* M, LAPACK_CONST blas_size_t* N, 
	blas_double_t* A, LAPACK_CONST blas_size_t* lda,
	blas_size_t* ipiv, blas_size_t* info);
extern "C" lapack_ret_t
cgetrf_(LAPACK_CONST blas_size_t* M, LAPACK_CONST blas_size_t* N, 
	blas_scomplex_t* A, LAPACK_CONST blas_size_t* lda,
	blas_size_t* ipiv, blas_size_t* info);
extern "C" lapack_ret_t
zgetrf_(LAPACK_CONST blas_size_t* M, LAPACK_CONST blas_size_t* N, 
	blas_complex_t* A, LAPACK_CONST blas_size_t* lda,
	blas_size_t* ipiv, blas_size_t* info);
    
// Getrs (solves system of equations from LU Decomposition)
extern "C" lapack_ret_t
sgetrs_(LAPACK_CONST char* trans, LAPACK_CONST blas_size_t* n, 
	LAPACK_CONST blas_size_t* n_rhs, LAPACK_CONST blas_float_t* A,
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_size_t* ipiv,
	blas_float_t* B, LAPACK_CONST blas_size_t* ldb,
	blas_size_t* info);
extern "C" lapack_ret_t
dgetrs_(LAPACK_CONST char* trans, LAPACK_CONST blas_size_t* n, 
	LAPACK_CONST blas_size_t* n_rhs, LAPACK_CONST blas_double_t* A,
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_size_t* ipiv,
	blas_double_t* B, LAPACK_CONST blas_size_t* ldb,
	blas_size_t* info);
extern "C" lapack_ret_t
cgetrs_(LAPACK_CONST char* trans, LAPACK_CONST blas_size_t* n, 
	LAPACK_CONST blas_size_t* n_rhs, LAPACK_CONST blas_scomplex_t* A,
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_size_t* ipiv,
	blas_scomplex_t* B, LAPACK_CONST blas_size_t* ldb,
	blas_size_t* info);
extern "C" lapack_ret_t
zgetrs_(LAPACK_CONST char* trans, LAPACK_CONST blas_size_t* n, 
	LAPACK_CONST blas_size_t* n_rhs, LAPACK_CONST blas_complex_t* A,
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_size_t* ipiv,
	blas_complex_t* B, LAPACK_CONST blas_size_t* ldb,
	blas_size_t* info);

// Getri (computes inverse from LU Decomposition)
extern "C" lapack_ret_t
sgetri_(LAPACK_CONST blas_size_t* N, blas_float_t* A, 
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_size_t* IPIV, 
	blas_float_t* WORK, LAPACK_CONST blas_size_t* lwork, 
	blas_size_t* INFO);
extern "C" lapack_ret_t
dgetri_(LAPACK_CONST blas_size_t* N, blas_double_t* A, 
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_size_t* IPIV, 
	blas_double_t* WORK, LAPACK_CONST blas_size_t* lwork, 
	blas_size_t* INFO);
extern "C" lapack_ret_t
cgetri_(LAPACK_CONST blas_size_t* N, blas_scomplex_t* A, 
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_size_t* IPIV, 
	blas_scomplex_t* WORK, LAPACK_CONST blas_size_t* lwork, 
	blas_size_t* INFO);
extern "C" lapack_ret_t
zgetri_(LAPACK_CONST blas_size_t* N, blas_complex_t* A, 
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_size_t* IPIV, 
	blas_complex_t* WORK, LAPACK_CONST blas_size_t* lwork, 
	blas_size_t* INFO);

//////////////////////////
// QR Decomposition
// Geqrf (does the QR Decomposition)
extern "C" lapack_ret_t
sgeqrf_(LAPACK_CONST blas_size_t* m, LAPACK_CONST blas_size_t* n,	
	blas_float_t* A, LAPACK_CONST blas_size_t* lda,
	blas_float_t* tau, blas_float_t* work, 
	LAPACK_CONST blas_size_t* lwork, blas_size_t* info);
extern "C" lapack_ret_t
dgeqrf_(LAPACK_CONST blas_size_t* m, LAPACK_CONST blas_size_t* n,	
	blas_double_t* A, LAPACK_CONST blas_size_t* lda,
	blas_double_t* tau, blas_double_t* work, 
	LAPACK_CONST blas_size_t* lwork, blas_size_t* info);
extern "C" lapack_ret_t
cgeqrf_(LAPACK_CONST blas_size_t* m, LAPACK_CONST blas_size_t* n,	
	blas_scomplex_t* A, LAPACK_CONST blas_size_t* lda,
	blas_scomplex_t* tau, blas_scomplex_t* work, 
	LAPACK_CONST blas_size_t* lwork, blas_size_t* info);
extern "C" lapack_ret_t
zgeqrf_(LAPACK_CONST blas_size_t* m, LAPACK_CONST blas_size_t* n,	
	blas_complex_t* A, LAPACK_CONST blas_size_t* lda,
	blas_complex_t* tau, blas_complex_t* work, 
	LAPACK_CONST blas_size_t* lwork, blas_size_t* info);

// Orgqr / Ungqr (retrieves the matrix Q of the QR Decomposition)
extern "C" lapack_ret_t
sorgqr_(LAPACK_CONST blas_size_t* m, LAPACK_CONST blas_size_t* n,
	LAPACK_CONST blas_size_t* k, blas_float_t* A, 
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_float_t* tau, 
	blas_float_t* work, LAPACK_CONST blas_size_t* lwork, 
	blas_size_t* info);
extern "C" lapack_ret_t
dorgqr_(LAPACK_CONST blas_size_t* m, LAPACK_CONST blas_size_t* n,
	LAPACK_CONST blas_size_t* k, blas_double_t* A, 
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_double_t* tau, 
	blas_double_t* work, LAPACK_CONST blas_size_t* lwork, 
	blas_size_t* info);
extern "C" lapack_ret_t
cungqr_(LAPACK_CONST blas_size_t* m, LAPACK_CONST blas_size_t* n,
	LAPACK_CONST blas_size_t* k, blas_scomplex_t* A, 
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_scomplex_t* tau, 
	blas_scomplex_t* work, LAPACK_CONST blas_size_t* lwork, 
	blas_size_t* info);
extern "C" lapack_ret_t
zungqr_(LAPACK_CONST blas_size_t* m, LAPACK_CONST blas_size_t* n,
	LAPACK_CONST blas_size_t* k, blas_complex_t* A, 
	LAPACK_CONST blas_size_t* lda, LAPACK_CONST blas_complex_t* tau, 
	blas_complex_t* work, LAPACK_CONST blas_size_t* lwork, 
	blas_size_t* info);


//////////////////////////
// Cholesky Decomposition
extern "C" lapack_ret_t
spotrf_(LAPACK_CONST char* uplo, LAPACK_CONST blas_size_t* n, blas_float_t* A,
	LAPACK_CONST blas_size_t* lda, blas_size_t* info);

extern "C" lapack_ret_t
dpotrf_(LAPACK_CONST char* uplo, LAPACK_CONST blas_size_t* n, blas_double_t* A,
	LAPACK_CONST blas_size_t* lda, blas_size_t* info);

extern "C" lapack_ret_t
cpotrf_(LAPACK_CONST char* uplo, LAPACK_CONST blas_size_t* n, blas_scomplex_t* A,
	LAPACK_CONST blas_size_t* lda, blas_size_t* info);

extern "C" lapack_ret_t
zpotrf_(LAPACK_CONST char* uplo, LAPACK_CONST blas_size_t* n, blas_complex_t* A,
	LAPACK_CONST blas_size_t* lda, blas_size_t* info);

//////////////////////////
// Eigenvalues

// symmetric/hermitian
extern "C" lapack_ret_t
ssyev_(LAPACK_CONST char* jobz, LAPACK_CONST char* uplo, 
       LAPACK_CONST blas_size_t* n, blas_float_t* a, 
       LAPACK_CONST blas_size_t* lda, blas_float_t* w, 
       blas_float_t* work, LAPACK_CONST blas_size_t* lwork, 
       blas_size_t* info);
    
extern "C" lapack_ret_t
dsyev_(LAPACK_CONST char* jobz, LAPACK_CONST char* uplo, 
       LAPACK_CONST blas_size_t* n, blas_double_t* a, 
       LAPACK_CONST blas_size_t* lda, blas_double_t* w, 
       blas_double_t* work, LAPACK_CONST blas_size_t* lwork, 
       blas_size_t* info);

extern "C" lapack_ret_t
cheev_(LAPACK_CONST char* jobz, LAPACK_CONST char* uplo, 
       LAPACK_CONST blas_size_t* n, blas_scomplex_t* a, 
       LAPACK_CONST blas_size_t* lda, blas_float_t* w, 
       blas_scomplex_t* work, LAPACK_CONST blas_size_t* lwork, 
       blas_float_t* rwork, blas_size_t* info);
    
extern "C" lapack_ret_t
zheev_(LAPACK_CONST char* jobz, LAPACK_CONST char* uplo, 
       LAPACK_CONST blas_size_t* n, blas_complex_t* a, 
       LAPACK_CONST blas_size_t* lda, blas_double_t* w, 
       blas_complex_t* work, LAPACK_CONST blas_size_t* lwork, 
       blas_double_t* rwork, blas_size_t* info);

// generic
extern "C" lapack_ret_t
sgeev_(LAPACK_CONST char* jobvl, LAPACK_CONST char* jobvr,
       LAPACK_CONST blas_size_t* n, blas_float_t* a, 
       LAPACK_CONST blas_size_t* lda, blas_float_t* wr, 
       blas_float_t* wi, blas_float_t* vl,
       LAPACK_CONST blas_size_t* ldvl, blas_float_t* vr,
       LAPACK_CONST blas_size_t* ldvr, blas_float_t* work,
       LAPACK_CONST blas_size_t* lwork, blas_size_t* info);
extern "C" lapack_ret_t
dgeev_(LAPACK_CONST char* jobvl, LAPACK_CONST char* jobvr,
       LAPACK_CONST blas_size_t* n, blas_double_t* a, 
       LAPACK_CONST blas_size_t* lda, blas_double_t* wr, 
       blas_double_t* wi, blas_double_t* vl,
       LAPACK_CONST blas_size_t* ldvl, blas_double_t* vr,
       LAPACK_CONST blas_size_t* ldvr, blas_double_t* work,
       LAPACK_CONST blas_size_t* lwork, blas_size_t* info);
extern "C" lapack_ret_t
cgeev_(LAPACK_CONST char* jobvl, LAPACK_CONST char* jobvr,
       LAPACK_CONST blas_size_t* n, blas_scomplex_t* a, 
       LAPACK_CONST blas_size_t* lda, blas_scomplex_t* w, 
       blas_scomplex_t* vl, LAPACK_CONST blas_size_t* ldvl, 
       blas_scomplex_t* vr, LAPACK_CONST blas_size_t* ldvr, 
       blas_scomplex_t* work, LAPACK_CONST blas_size_t* lwork, 
       blas_float_t* rwork, blas_size_t* info);
extern "C" lapack_ret_t
zgeev_(LAPACK_CONST char* jobvl, LAPACK_CONST char* jobvr,
       LAPACK_CONST blas_size_t* n, blas_complex_t* a, 
       LAPACK_CONST blas_size_t* lda, blas_complex_t* w, 
       blas_complex_t* vl, LAPACK_CONST blas_size_t* ldvl, 
       blas_complex_t* vr, LAPACK_CONST blas_size_t* ldvr, 
       blas_complex_t* work, LAPACK_CONST blas_size_t* lwork, 
       blas_double_t* rwork, blas_size_t* info);
    
// Real symmetric tridiagonal eigensolvers
extern "C" lapack_ret_t
sstev_(LAPACK_CONST char* jobz, LAPACK_CONST blas_size_t* N, 
       blas_float_t* D, blas_float_t* E, blas_float_t* Z,
       LAPACK_CONST blas_size_t* ldz, blas_float_t* work, 
       blas_size_t* info);

extern "C" lapack_ret_t
dstev_(LAPACK_CONST char* jobz, LAPACK_CONST blas_size_t* N, 
       blas_double_t* D, blas_double_t* E, blas_double_t* Z,
       LAPACK_CONST blas_size_t* ldz, blas_double_t* work, 
       blas_size_t* info);


///////////////////////////////////
// Generalized eigenvalue problems
extern "C" lapack_ret_t
ssygv_(LAPACK_CONST blas_size_t* itype, LAPACK_CONST char* jobz, 
       LAPACK_CONST char* uplo, LAPACK_CONST blas_size_t* N, 
       blas_float_t* A, LAPACK_CONST blas_size_t* lda, 
       blas_float_t* B, LAPACK_CONST blas_size_t* ldb,
       blas_float_t* W, blas_float_t* work, 
       LAPACK_CONST blas_size_t* lwork, blas_size_t* info);

extern "C" lapack_ret_t
dsygv_(LAPACK_CONST blas_size_t* itype, LAPACK_CONST char* jobz, 
       LAPACK_CONST char* uplo, LAPACK_CONST blas_size_t* N, 
       blas_double_t* A, LAPACK_CONST blas_size_t* lda, 
       blas_double_t* B, LAPACK_CONST blas_size_t* ldb, 
       blas_double_t* W, blas_double_t* work, 
       LAPACK_CONST blas_size_t* lwork, blas_size_t* info);

extern "C" lapack_ret_t
chegv_(LAPACK_CONST blas_size_t* itype, LAPACK_CONST char* jobz, 
       LAPACK_CONST char* uplo, LAPACK_CONST blas_size_t* N, 
       blas_scomplex_t* A, LAPACK_CONST blas_size_t* lda, 
       blas_scomplex_t* B, LAPACK_CONST blas_size_t* ldb,
       blas_float_t* W, blas_scomplex_t* work, 
       LAPACK_CONST blas_size_t* lwork, blas_float_t* rwork, 
       blas_size_t* info);
    
extern "C" lapack_ret_t
zhegv_(LAPACK_CONST blas_size_t* itype, LAPACK_CONST char* jobz, 
       LAPACK_CONST char* uplo, LAPACK_CONST blas_size_t* N, 
       blas_complex_t* A, LAPACK_CONST blas_size_t* lda, 
       blas_complex_t* B, LAPACK_CONST blas_size_t* ldb, 
       blas_double_t* W, blas_complex_t* work, 
       LAPACK_CONST blas_size_t* lwork, blas_double_t* rwork, 
       blas_size_t* info);

#endif

#endif
