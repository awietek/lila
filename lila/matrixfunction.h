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

#ifndef LILA_MATRIX_FUNCTION_H_
#define LILA_MATRIX_FUNCTION_H_

#include <math.h>

#include "matrix.h"
#include "mult.h"
#include "vector.h"
#include "eigen.h"
#include "complex.h"
#include "blaslapack.h"

namespace lila {

  template<class coeff_t, class function_t>
  inline void FunctionH(Matrix<coeff_t>& matrix, function_t fun, char uplo='U')
  {
    using size_type = blaslapack::blas_size_t;

    assert(matrix.nrows() == matrix.ncols());
    const size_type n = matrix.nrows();
    Matrix<coeff_t> Qmat = matrix;
    Vector<real_t<coeff_t>> eigs = EigenHDestroy(Qmat, true, uplo);
    Matrix<coeff_t> Bmat = Qmat;

    // Multiply with diagonal matrix, where function has been applied
    int incx = 1;    
    for(int j : eigs.rows())
      {
	coeff_t fun_of_eig = static_cast<coeff_t>(eigs(j));
	fun(fun_of_eig);

	// Scale j-th column
	blaslapack::scal(&n, &fun_of_eig, Bmat.data() + j*n, &incx); 
      }
    Zeros(matrix);
    Mult(Bmat, Qmat, matrix, (coeff_t)1., (coeff_t)0., 'N', 'C');
  }
  template<class coeff_t>
  inline void ExpH(Matrix<coeff_t>& matrix, coeff_t alpha=1., char uplo='U')
  { 
    FunctionH(matrix, 
	      [alpha](coeff_t& x) { 
		x = exp(alpha * x);
	      }, 
	      uplo); 
  }
  
  template<class coeff_t>
  inline void LogH(Matrix<coeff_t>& matrix, coeff_t alpha=1., char uplo='U')
  { 
    FunctionH(matrix, 
	      [alpha](coeff_t& x) { 
		x = log(x) / alpha; 
	      },
	      uplo); 
  }

}

#endif
