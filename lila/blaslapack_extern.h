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

#ifndef PLATFORM_MKL

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
    extern "C" blas_scomplex_t cdotu_(const blas_size_t* N,
				      const blas_scomplex_t* x,
				      const blas_size_t* incx, 
				      const blas_scomplex_t* y,
				      const blas_size_t* incy);
    extern "C" blas_complex_t zdotu_(const blas_size_t* N,
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
    extern "C" void zgemv_(const char* trans, const blas_size_t* m,c
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
    

  }  // namespace lila
}  // namespace blaslapack

#endif

#endif
