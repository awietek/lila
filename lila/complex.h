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

#ifndef LILA_COMPLEX_H_
#define LILA_COMPLEX_H_

#include <complex>

#include "matrix.h"
#include "vector.h"
#include "detail/complex_detail.h"

namespace lila {
  
  template <class coeff_t>
  using real_t = typename detail::real_type_struct<coeff_t>::type;

  template <class coeff_t>
  using complex_t = typename detail::complex_type_struct<coeff_t>::type;

  inline float real (float x) { return x; }
  inline double real (double x) { return x; }
  inline float real (std::complex<float> x) { return x.real(); }
  inline double real (std::complex<double> x) { return x.real(); }

  inline float imag (float) { return 0.; }
  inline double imag (double) { return 0.; }
  inline float imag (std::complex<float> x) { return x.imag(); }
  inline double imag (std::complex<double> x) { return x.imag(); }

  inline float conj (float x) { return x; }
  inline double conj (double x) { return x; }
  inline std::complex<float> conj (std::complex<float> x) { return std::conj(x); }
  inline std::complex<double> conj (std::complex<double> x) { return std::conj(x); }

  template <class coeff_t>
  inline Matrix<real_t<coeff_t>> Real(const Matrix<coeff_t>& X)
  {
    Matrix<real_t<coeff_t>> Y(X.nrows(), Y.ncols());
    for (auto i : X.rows())
      for (auto j : X.cols())
	Y(i, j) = real(X(i, j));
    return Y;
  }

  template <class coeff_t>
  inline Vector<real_t<coeff_t>> Real(const Vector<coeff_t>& X)
  {
    Vector<real_t<coeff_t>> Y(X.nrows());
    for (auto i : X.rows())
	Y(i) = real(X(i));
    return Y;
  }
 
  template <class coeff_t>
  inline Matrix<real_t<coeff_t>> Imag(const Matrix<coeff_t>& X)
  {
    Matrix<real_t<coeff_t>> Y(X.nrows(), Y.ncols());
    for (auto i : X.rows())
      for (auto j : X.cols())
	Y(i, j) = imag(X(i, j));
    return Y;
  }

  template <class coeff_t>
  inline Vector<real_t<coeff_t>> Imag(const Vector<coeff_t>& X)
  {
    Vector<real_t<coeff_t>> Y(X.nrows());
    for (auto i : X.rows())
	Y(i) = imag(X(i));
    return Y;
  }  
  
  template <class coeff_t>
  Matrix<coeff_t> Conj(const Matrix<coeff_t>& mat) 
  {
    Matrix<coeff_t> mat_c(mat.nrows(), mat.ncols());
    for (auto i : mat.rows())
      for (auto j : mat.cols())
	mat_c(i, j) = conj(mat(i, j));
    return mat_c;
  }

  template <class coeff_t>
  Matrix<coeff_t> Herm(const Matrix<coeff_t>& mat) 
  {
    return Conj(Transpose(mat));
  }

}

#endif  // LILA_COMPLEX_H_
