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

#include "blaslapack_types.h"

#ifdef PLATFORM_MKL  // Using the Intel MKL

#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include "mkl_types.h"
#include "mkl.h"
#define __LAPACK_ROUTINE_NAME(x) x

#else  // Using normal LAPACK  
#include "blaslapack_extern.h"
#define __LAPACK_ROUTINE_NAME(x) x##_
#endif

namespace lila { 
  namespace blaslapack {

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
    { __LAPACK_ROUTINE_NAME(sdot)(N, x, incx, y, incy); }
    inline double dot(const blas_size_t* N ,const blas_double_t* x,
		      const blas_size_t* incx, const blas_double_t* y,
		      const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(ddot)(N, x, incx, y, incy); }
    inline blas_scomplex_t dot(const blas_size_t* N, const blas_scomplex_t* x,
			      const blas_size_t* incx, const blas_scomplex_t* y,
			      const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(cdotu)(N, x, incx, y, incy); }
    inline blas_complex_t dot(const blas_size_t* N, const blas_complex_t* x,
			      const blas_size_t* incx, const blas_complex_t* y,
			      const blas_size_t* incy)
    { __LAPACK_ROUTINE_NAME(zdotu)(N, x, incx, y, incy); }
        
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
    inline void gemv(const char* trans, const blas_size_t* m,c
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

  }  // namespace lila
}  // namespace blaslapack

#endif
