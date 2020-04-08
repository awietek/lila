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

#ifndef LILA_MULT_H_
#define LILA_MULT_H_

#include "matrix.h"
#include "vector.h"
#include "blaslapack/blaslapack.h"
#include "special.h" 
#include "complex.h"
#include "range.h"

namespace lila {

  template <class coeff_t>
  inline void Mult(const Matrix<coeff_t>& A, const Matrix<coeff_t>& B, 
		   Matrix<coeff_t>& C, coeff_t alpha = 1., coeff_t beta = 0., 
		   char transa = 'N', char transb = 'N')
  {
    using size_type = blaslapack::blas_size_t;

    // matrix dimensions
    size_type m = (transa == 'N') ? A.nrows() : A.ncols();
    size_type n = (transb == 'N') ? B.ncols() : B.nrows();
    size_type ka = (transa == 'N') ? A.ncols() : A.nrows();
    size_type kb = (transb == 'N') ? B.nrows() : B.ncols();

    assert(ka == kb); // Check if valid multiplication dimensions

    if ((C.nrows() != m) || (C.ncols() != n)) C.resize(m, n);

    // leading dimensions
    size_type lda = (transa == 'N') ? m : ka;
    size_type ldb = (transb == 'N') ? ka : n;
    size_type ldc = m;

    blaslapack::gemm(&transa,
		     &transb,
		     &m,
		     &n,
		     &ka,
		     LILA_BLAS_CAST(coeff_t,&alpha),
		     LILA_BLAS_CONST_CAST(coeff_t,A.data()),
		     &lda,
		     LILA_BLAS_CONST_CAST(coeff_t,B.data()),
		     &ldb,
		     LILA_BLAS_CAST(coeff_t,&beta),
		     LILA_BLAS_CAST(coeff_t,C.data()),
		     &ldc);
  }
  
  
  template <class coeff_t>
  inline Matrix<coeff_t> Mult(const Matrix<coeff_t>& A, const Matrix<coeff_t>& B)
  {
    Matrix<coeff_t> C;
    Mult(A, B, C);
    return C;
  }

  template <class coeff_t>
  inline void Mult(const Matrix<coeff_t>& A, const Vector<coeff_t>& X, 
		   Vector<coeff_t>& Y, coeff_t alpha = 1., coeff_t beta = 0., 
		   char trans = 'N')
  {
    using size_type = blaslapack::blas_size_t;

    // matrix dimensions
    size_type m = (trans == 'N') ? A.nrows() : A.ncols();
    size_type n = (trans == 'N') ? A.ncols() : A.nrows();

    assert(n == X.size()); // Check if valid multiplication dimensions

    if (Y.size() != m) Y.resize(n);

    // leading dimensions
    size_type lda = A.nrows();
    size_type inc = 1;
    blaslapack::gemv(&trans,
		     &m,
		     &n,
		     LILA_BLAS_CAST(coeff_t,&alpha),
		     LILA_BLAS_CONST_CAST(coeff_t,A.data()),
		     &lda,
		     LILA_BLAS_CONST_CAST(coeff_t,X.data()),
		     &inc, 
		     LILA_BLAS_CAST(coeff_t,&beta),
		     LILA_BLAS_CAST(coeff_t,Y.data()),
		     &inc);
  }
  
  template <class coeff_t>
  inline Vector<coeff_t> Mult(const Matrix<coeff_t>& A, const Vector<coeff_t>& X)
  {
    Vector<coeff_t> Y;
    Mult(A, X, Y);
    return Y;
  }

  template <class coeff_t>
  inline void Kron(const Matrix<coeff_t>& A, const Matrix<coeff_t>& B, 
		   Matrix<coeff_t>& C)
  {
    using size_type = blaslapack::blas_size_t;
    size_type m = A.nrows();
    size_type n = A.ncols();
    size_type p = B.nrows();
    size_type q = B.ncols();

    size_type mc = m*p;
    size_type nc = n*q;
    C.resize(mc, nc);
    for (int s : range<int>(n))
      for (int r : range<int>(m))
	for (int w : range<int>(q))
	  for (int v : range<int>(p))
	    C(p*r+v, q*s+w) = A(r, s)*B(v,w);
  }


}

#endif
