// Copyright 2019 Alexander Wietek - All Rights Reserved.
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

#ifndef LILA_CHOLESKY_H_
#define LILA_CHOLESKY_H_

#include "matrix.h"
#include "blaslapack/blaslapack.h"

namespace lila {

  /*! @brief Computes Cholesky decomposition inplace (destroying input)

    Cholesky decomposition of a positive definite symmetric (hermitian) matrix A
    \f$ A =  U**H * U \f$ if uplo = 'U'
    \f$ A =  L * L**H \f$ if uplo = 'L'
    where U is upper triangular, and L lower triangular

    @param A lila::Matrix (overwritten, will contain results of potrf LAPACK call 
                           (warning: this is in general neither U nor L!)) 
    @tparam coeff_t type of coefficients of object 
   */
  template <class coeff_t>
  inline void 
  CholeskyInplace(Matrix<coeff_t>& A, char uplo = 'U')
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == A.ncols());
    const size_type n = A.nrows();
    const size_type lda = n;
    int info = 0;
    blaslapack::potrf(&uplo, &n, A.data(), &lda, &info);
  }

  /*! @brief Computes Cholesky decomposition

    Cholesky decomposition of a positive definite symmetric (hermitian) matrix A
    \f$ A =  U**H * U \f$ if uplo = 'U'
    \f$ A =  L * L**H \f$ if uplo = 'L'
    where U is upper triangular, and L lower triangular

    @param A lila::Matrix (overwritten, contains factors U or L respectively) 
    @tparam coeff_t type of coefficients of object 

    @return U/L the upper/lower triangular matrix of the Cholesky decomposition
   */
  template <class coeff_t>
  inline Matrix<coeff_t> 
  Cholesky(const Matrix<coeff_t>& A, char uplo = 'U')
  {
    auto res = A;
    CholeskyInplace(res, uplo);

    // Set lower triangular to zero
    int n = A.nrows();
    if (uplo == 'U')
      {
	for (int j=0; j<n; ++j)	
	  for (int i=j+1; i<n; ++i)
	    res(i, j) = 0;
      }
    else if (uplo == 'L')
      {
	for (int i=0; i<n; ++i)	
	  for (int j=i+1; j<n; ++j)
	    res(i, j) = 0;
      }
    return res;
  }


  

}

#endif
