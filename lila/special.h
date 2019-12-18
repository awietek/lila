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
#include "vector.h"
#include "range.h"
#include "complex.h"
#include "mult.h"
#include "compare.h"
#include "print.h"

namespace lila {

  template <class coeff_t>
  inline void Zeros(Vector<coeff_t>& vec)
  { std::fill(vec.begin(), vec.end(), 0.); }

  template <class coeff_t>
  inline void Zeros(Matrix<coeff_t>& mat)
  { std::fill(mat.begin(), mat.end(), 0.); }

  template <class coeff_t>
  Vector<coeff_t> Zeros(size_type m)
  { 
    Vector<coeff_t> vec(m);
    Zeros(vec);
    return vec;
  }

  template <class coeff_t>
  Vector<coeff_t> ZerosLike(const Vector<coeff_t>& vec)
  { 
    Vector<coeff_t> res(vec.size());
    Zeros(res);
    return res;
  }

  template <class coeff_t>
  Matrix<coeff_t> Zeros(size_type m, size_type n)
  { 
    Matrix<coeff_t> mat(m, n);
    Zeros(mat);
    return mat;
  }

  template <class coeff_t>
  Matrix<coeff_t> ZerosLike(const Matrix<coeff_t>& mat)
  { 
    Matrix<coeff_t> res(mat.nrows(), mat.ncols());
    Zeros(res);
    return res;
  }

  template <class matrix_t>
  inline void Ones(matrix_t& mat)
  { std::fill(mat.begin(), mat.end(), 1.); }
  
  template <class coeff_t>
  inline void Identity(Matrix<coeff_t>& mat)
  { 
    Zeros(mat);
    for (auto i : range<size_type>(std::min(mat.nrows(), mat.ncols())))
      mat(i,i) = 1.;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> Identity(size_type m)
  { 
    auto id = Zeros<coeff_t>(m, m);
    Identity(id);
    return id;
  }

  template <class coeff_t>
  inline Vector<complex_t<coeff_t>> Complex(Vector<coeff_t>& vec)
  { 
    auto complex_vec = Zeros<complex_t<coeff_t>>(vec.nrows());
    for (auto j : vec.rows())
      complex_vec(j) = (complex_t<coeff_t>)vec(j);
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
  Vector<coeff_t> linspace(coeff_t a, coeff_t b, int n) {
    Vector<coeff_t> vector(n);
    coeff_t step = (b-a) / (n-1);
    for (int i=0; i<n; ++i)
      vector(i) = a + step*i;
    return vector;
  }
}

#endif
