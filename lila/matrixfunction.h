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

#include <lila/add.h>
#include <lila/matrix.h>
#include <lila/mult.h>
#include <lila/vector.h>
#include <lila/eigen.h>
#include <lila/complex.h>
#include <lila/blaslapack/blaslapack.h>
#include <lila/special.h>
#include <lila/solve.h>

namespace lila {

  template<class coeff_t, class function_t>
  inline void FunctionSym(Matrix<coeff_t>& matrix, function_t fun, char uplo='U')
  {
    using size_type = blaslapack::blas_size_t;

    assert(matrix.nrows() == matrix.ncols());
    const size_type n = matrix.nrows();
    Matrix<coeff_t> Qmat = matrix;
    Vector<real_t<coeff_t>> eigs = EigenSymInplace(Qmat, true, uplo);
    Matrix<coeff_t> Bmat = Qmat;

    // Multiply with diagonal matrix, where function has been applied
    int incx = 1;    
    for(int j : eigs.rows())
      {
	coeff_t fun_of_eig = static_cast<coeff_t>(eigs(j));
	fun(fun_of_eig);

	// Scale j-th column
	blaslapack::scal(&n,
			 LILA_BLAS_CAST(coeff_t,&fun_of_eig),
			 LILA_BLAS_CAST(coeff_t,Bmat.data()) + j*n,
			 &incx); 
      }
    Zeros(matrix);
    Mult(Bmat, Qmat, matrix, (coeff_t)1., (coeff_t)0., 'N', 'C');
  }

  template<class coeff_t>
  inline void ExpSym(Matrix<coeff_t>& matrix, coeff_t alpha=1., char uplo='U')
  { 
    FunctionSym(matrix, 
		[alpha](coeff_t& x) { 
		  x = exp(alpha * x);
		}, 
		uplo); 
  }
  
  template<class coeff_t>
  inline void LogSym(Matrix<coeff_t>& matrix, coeff_t alpha=1., char uplo='U')
  { 
    FunctionSym(matrix, 
		[alpha](coeff_t& x) { 
		  x = log(x) / alpha; 
		},
		uplo); 
  }
  
  /*! @brief Computes exponential of generic matrix (EXPOKIT routine)
    
    the exponential exp(alpha*A) of a matrix, where alpha is a prefactor
    
    @param A lila::Matrix
    @param alpha prefector in the exponent (optional: default is 1.)
    @tparam coeff_t type of coefficients of matrix 
    
    @return exponential exp(alpha*A) of matrix
  */
  template <class coeff_t>
  inline Matrix<coeff_t> ExpM(const Matrix<coeff_t>& A, coeff_t alpha=1.)
  {
    int n = A.nrows();
    assert(n == A.ncols());

    const int q = 6;
    Matrix<coeff_t> a2 = alpha * A;
    real_t<coeff_t> a_norm = NormLi(a2);
    int ee = (int)log2abs(a_norm) + 1;
    int s = std::max(0, ee + 1);
    coeff_t t = 1.0 / pow ( 2.0, s );
    Scale(t, a2);
    Matrix<coeff_t> x = a2;
    coeff_t c = 0.5;
    Matrix<coeff_t> e = Identity<coeff_t>(n);
    Add(a2, e, c);
    Matrix<coeff_t> d = Identity<coeff_t>(n);
    Add(a2, d, -c);
    int p = 1;

    for (int k = 2; k <= q; k++ )
      {
	c = c * (coeff_t) ( q - k + 1 ) / (coeff_t) ( k * ( 2 * q - k + 1 ) );

	x = Mult(a2, x);
	Add(x, e, c);
	if ( p )
	  {
	    Add(x, d, c);
	  }
	else
	  {
	    Add(x, d, -c);
	  }
	p = !p;
      }
    Solve(d, e);
    for (int k = 1; k <= s; k++ )
      {
	e = Mult(e, e);
      }
    return e;
  }



  template <class coeff_t>
  inline Matrix<coeff_t> Unitary(size_type n, 
				 const Vector<real_t<coeff_t>>& params)
  { 
    assert(n>0);
    assert(params.size() == n*n);
    
    Matrix<coeff_t> unitary(n, n);
    Zeros(unitary);
    
    // Set diagonal elements
    for (int i : range<int>(n))
      unitary(i, i) = params(i);
    
    int param_offset = 0;
    for (int n_offdiag=1; n_offdiag < n; ++n_offdiag)
      for (int j : range<int>(n - n_offdiag))
    	{	  
    	  // Construct hermitian elementary matrix
    	  coeff_t entry = coeff_t(params(n + param_offset), 
    				  params(n + param_offset + 1));
	  
    	  unitary(j, n_offdiag + j) = entry;
    	  unitary(n_offdiag + j, j) = lila::conj(entry);

    	  // // IMPROVEMENT: do exponentiation without full matrix exponentiation
	  // Zeros(elementary_2x2);
	  // elementary_2x2(0, 1) = entry;
	  // elementary_2x2(0, 1) = lila::conj(entry);
	  // ExpH(elementary_2x2, coeff_t(0., 1.), 'U');

    	  // elementary(j, j) = elementary_2x2(0, 0);
    	  // elementary(n_offdiag + j, n_offdiag + j) = elementary_2x2(1, 1);
	  // elementary(j, n_offdiag + j) = elementary_2x2(0, 1);
    	  // elementary(n_offdiag + j, j) = elementary_2x2(1, 0);

    	  // Matrix<coeff_t> unitary_copy = unitary;
    	  // Mult(unitary_copy, elementary, unitary);
    	  param_offset += 2;
    	}   
    ExpSym(unitary, coeff_t(0., 1.), 'U');

    assert(param_offset == n*(n-1));
    return unitary;
  }

  
  
  

  template <class coeff_t>
  inline bool IsUnitary(const Matrix<coeff_t>& mat)
  {
    Matrix<coeff_t> p(mat);
    Zeros(p);
    Mult(mat, mat, p, coeff_t(1.), coeff_t(0.), 'C', 'N');
    Matrix<coeff_t> Id(p);
    Identity(Id);
    return close(Id, p);
  }

}

#endif
