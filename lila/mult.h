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
#include "blaslapack.h"
#include "special.h" 

namespace lila {

  template <class coeff_t>
  inline void Mult(const Matrix<coeff_t>& A, const Matrix<coeff_t>& B, 
		   Matrix<coeff_t>& C, coeff_t alpha = 1., coeff_t beta = 0., 
		   char transa = 'N', char transb = 'N')
  {
    using size_type = blaslapack::blas_size_t;

    // matrix dimensions
    const size_type m = (transa == 'N') ? A.nrows() : A.ncols();
    const size_type n = (transb == 'N') ? B.ncols() : B.nrows();
    const size_type ka = (transa == 'N') ? A.ncols() : A.nrows();
    const size_type kb = (transb == 'N') ? B.nrows() : B.ncols();

    assert(ka == kb); // Check if valid multiplication dimensions

    if ((C.nrows() != m) || (C.ncols() != n)) C.resize(m, n);

    // leading dimensions
    const size_type lda = (transa == 'N') ? m : ka;
    const size_type ldb = (transb == 'N') ? ka : n;
    const size_type ldc = m;

    blaslapack::gemm(&transa, &transb, &m, &n, &ka, &alpha,
		     A.data(), &lda, B.data(), &ldb, &beta, C.data(), &ldc);
  }

  template <class coeff_t>
  inline void Mult(const Matrix<coeff_t>& A, const Vector<coeff_t>& X, 
		   Vector<coeff_t>& Y, coeff_t alpha = 1., coeff_t beta = 0., 
		   char trans = 'N')
  {
    using size_type = blaslapack::blas_size_t;

    // matrix dimensions
    const size_type m = (trans == 'N') ? A.nrows() : A.ncols();
    const size_type n = (trans == 'N') ? A.ncols() : A.nrows();

    assert(n == X.size()); // Check if valid multiplication dimensions

    if (Y.size() != n) Y.resize(n);

    // leading dimensions
    const size_type lda = A.nrows();

    blaslapack::gemv(&transa, &transb, &m, &n, &ka, &alpha,
		     A.data(), &lda, B.data(), &ldb, &beta, C.data(), &ldc);
  }

}

#endif
