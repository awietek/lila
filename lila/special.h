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

#ifndef LILA_SPECIAL_H_
#define LILA_SPECIAL_H_

#include <algorithm>

#include "common.h"
#include "matrix.h"
#include "matrixfunction.h"
#include "vector.h"
#include "range.h"
#include "complex.h"
#include "mult.h"
#include "compare.h"
#include "print.h"

namespace lila {

  template <class matrix_t>
  inline void Zeros(matrix_t& mat)
  { std::fill(mat.begin(), mat.end(), 0.); }

  template <class coeff_t>
  Vector<coeff_t> Zeros(int m)
  { 
    Vector<coeff_t> vec(m);
    Zeros(vec);
    return vec;
  }

  template <class coeff_t>
  Matrix<coeff_t> Zeros(int m, int n)
  { 
    Matrix<coeff_t> mat(m, n);
    Zeros(mat);
    return mat;
  }

  template <class matrix_t>
  inline void Ones(matrix_t& mat)
  { std::fill(mat.begin(), mat.end(), 1.); }
  
  template <class matrix_t>
  inline void Identity(matrix_t& mat)
  { 
    Zeros(mat);
    for (auto i : range<int>(std::min(mat.nrows(), mat.ncols())))
      mat(i,i) = 1.;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> Identity(int m)
  { 
    auto id = Zeros<coeff_t>(m, m);
    Identity(id);
    return id;
  }

  template <class coeff_t>
  inline Vector<complex_t<coeff_t>> Complex(Vector<coeff_t>& vec)
  { 
    printf("tytyt %d\n", vec.nrows());
    auto complex_vec = Zeros<complex_t<coeff_t>>(vec.nrows());
    for (auto j : vec.rows())
      complex_vec(j) = (complex_t<coeff_t>)vec(j);
    printf("uiuiu\n");
    return complex_vec;
  }

  template <class coeff_t>
  inline Matrix<complex_t<coeff_t>> Complex(Matrix<coeff_t>& mat)
  { 
    auto complex_mat = Zeros<complex_t<coeff_t>>(mat.nrows(), mat.ncols());
    for (auto i : mat.rows())
      for (auto j : mat.cols())
	complex_mat(i, j) = (complex_t<coeff_t>)mat(i, j);
    return complex_mat;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> Diag(const Vector<coeff_t>& diag, int offset=0)
  { 
    int absoffset = std::abs(offset);
    int dim = diag.size() + std::abs(offset);
    Matrix<coeff_t> mat = Zeros<coeff_t>(dim, dim);
    for (int i : diag.rows())
      {
	if (offset < 0) mat(i + absoffset, i) = diag(i);
	else mat(i, i + absoffset) = diag(i);
      }
    return mat;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> Unitary(int n, 
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
    ExpH(unitary, coeff_t(0., 1.), 'U');

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


  template <class coeff_t>
  Vector<coeff_t> linspace(coeff_t a, coeff_t b, int n) {
    Vector<coeff_t> vector(n);
    coeff_t step = (b-a) / (n-1);
    for (int i=0; i<n; ++i)
      vector(i) = a + step*i;
    return vector;
  }
}

#endif
