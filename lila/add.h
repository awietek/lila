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

/*!
  @file add.h
  @brief Basic routines for manipulating Vectors and Matrices.
 */

#ifndef LILA_ADD_H_
#define LILA_ADD_H_

#include <cassert>
#include <algorithm>
#include <numeric>

#include "complex.h"
#include "matrix.h"
#include "vector.h"
#include "blaslapack.h"
#include "compare.h"

namespace lila {

  /*! @brief Maps a function onto a lila::Matrix or lila::Vector.

    __Usage example (to square entries of vector)__
    @code
      auto sq_vec = vec;
      lila::Map(sq_vec, [](double& e) { e = e*e; });
    @endcode

    @param X lila::Matrix / lila::Vector on which the function is applied
    @param func function object to be applied to the Matrix/Vector
    @tparam object_t either lila::Matrix or lila::Vector
    @tparam function_t function object type
    @return Matrix/Vector after function has been applied
   */
  template <class object_t, class function_t>
  inline object_t Map(object_t&& X, function_t func)
  { 
    std::for_each(X.begin(), X.end(), func);
    return X;
  }

  /*! @brief copies one lila::Matrix or lila::Vector to another

    @param X lila::Matrix / lila::Vector which is copied
    @param Y lila::Matrix / lila::Vector that is overwritten
    @tparam object_t either lila::Matrix or lila::Vector
    @tparam coeff_t type of coefficients of object 
   */
  template <class coeff_t>
  inline void Copy(const Vector<coeff_t>& X,  Vector<coeff_t>& Y)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::copy(&dx, X.data(), &inc, Y.data(), &inc);
  }

  template <class coeff_t>
  inline void Copy(const Matrix<coeff_t>& X,  Matrix<coeff_t>& Y)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::copy(&dx, X.data(), &inc, Y.data(), &inc);
  }

  /*! @brief Add a scalar multiple of lila::Matrix or lila::Vector to another

    Performs the following operation:
    \f$ \alpha X + Y \rightarrow Y \f$

    @param X lila::Matrix / lila::Vector 
    @param Y lila::Matrix / lila::Vector 
    @param alpha coeff_t scalar multiple of X

    @param func function object to be applied to the Matrix/Vector
    @tparam object_t either lila::Matrix or lila::Vector
    @tparam coeff_t type of coefficients of object 
   */
  template <class coeff_t>
  inline void Add(const Vector<coeff_t>& X, Vector<coeff_t>& Y, 
		  coeff_t alpha = static_cast<coeff_t>(1.))
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::axpy(&dx, &alpha, X.data(), &inc, Y.data(), &inc);
  }

  template <class coeff_t>
  inline void Add(const Matrix<coeff_t>& X,  Matrix<coeff_t>& Y, 
		  coeff_t alpha = static_cast<coeff_t>(1.))
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::axpy(&dx, &alpha, X.data(), &inc, Y.data(), &inc);
  }

  /*! @brief Multiplies object by a scalar

    Performs the following operation:
    \f$ \alpha X \rightarrow X \f$

    @param alpha coeff_t scalar multiple of X
    @param X lila::Matrix / lila::Vector which is scaled

    @param func function object to be applied to the Matrix/Vector
    @tparam object_t either lila::Matrix or lila::Vector
    @tparam coeff_t type of coefficients of object 
   */
  template <class coeff_t>
  inline void Scale(const coeff_t& alpha,  Vector<coeff_t>& X)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type inc = 1;
    blaslapack::scal(&dx, &alpha, X.data(), &inc);
  }
  template <class coeff_t>
  inline void Scale(const coeff_t& alpha,  Matrix<coeff_t>& X)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type inc = 1;
    blaslapack::scal(&dx, &alpha, X.data(), &inc);
  }


  /*! @brief Dot product between two objects
    
    Performs the following operation:
    \f$ \left< X | Y \right> \f$
    the vector X is conjugated in the complex case.

    @param X lila::Vector 
    @param Y lila::Vector
    @return dot product
    @tparam object_t either lila::Matrix or lila::Vector
    @tparam coeff_t type of coefficients of object 
   */
  template <class coeff_t>
  inline coeff_t Dot(const Vector<coeff_t>& X, const Vector<coeff_t>& Y)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    return blaslapack::dot(&dx, X.data(), &inc, Y.data(), &inc);
  }

  template <class coeff_t>
  inline coeff_t Dot(const Matrix<coeff_t>& X, const Matrix<coeff_t>& Y)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    return blaslapack::dot(&dx, X.data(), &inc, Y.data(), &inc);
  }



  template <class coeff_t>
  inline real_t<coeff_t> Norm(const Vector<coeff_t>& X)
  {
    // TODO: use proper LAPACK function here
    return sqrt(real(Dot(X,X)));
  }

  template <class coeff_t>
  inline real_t<coeff_t> Norm(const Matrix<coeff_t>& X)
  {
    // TODO: use proper LAPACK function here
    return sqrt(real(Dot(X,X)));
  }

  template <class coeff_t>
  inline coeff_t Sum(const Vector<coeff_t>& X)
  { return std::accumulate(X.begin(), X.end(), 0.); }
  template <class coeff_t>
  inline coeff_t Sum(const Matrix<coeff_t>& X)
  { return std::accumulate(X.begin(), X.end(), 0.); }


  
  template <class coeff_t>
  inline Vector<coeff_t> operator+
  (const Vector<coeff_t>& X, const Vector<coeff_t>& Y)
  {
    Vector<coeff_t> res(Y);
    Add(X, res);
    return res;
  }

  template <class coeff_t>
  inline Vector<coeff_t> operator-
  (const Vector<coeff_t>& X, const Vector<coeff_t>& Y)
  {
    Vector<coeff_t> res(X);
    Add(Y, res, static_cast<coeff_t>(-1.));
    return res;
  }

  template <class coeff_t>
  inline Vector<coeff_t> operator+
  (const Vector<coeff_t>& X, const coeff_t& c)
  {
    Vector<coeff_t> res(X);
    Map(res, [&c](coeff_t& x) { x = x + c; } );
    return res;
  }

  template <class coeff_t>
  inline Vector<coeff_t>& operator +=
  (Vector<coeff_t>& a, const Vector<coeff_t>& b)
  {
    a = a + b;
    return a;
  }
  template <class coeff_t>
  inline Vector<coeff_t>& operator -=
  (Vector<coeff_t>& a, const Vector<coeff_t>& b)
  {
    a = a - b;
    return a;
  }

  template <class coeff_t>
  inline Vector<coeff_t> operator-
  (const Vector<coeff_t>& X, const coeff_t& c)
  {
    Vector<coeff_t> res(X);
    Map(res, [&c](coeff_t& x) { x = x - c; } );
    return res;
  }


  template <class coeff_t>
  inline Vector<coeff_t> operator-(const Vector<coeff_t>& X)
  {
    Vector<coeff_t> res(X);
    Scale(static_cast<coeff_t>(-1.), res);
    return res;
  }

  template <class coeff_t>
  inline Vector<coeff_t> operator*
  (const coeff_t& alpha,  const Vector<coeff_t>& X)
  {
    Vector<coeff_t> res(X);
    Scale(alpha, res);
    return res;
  }

  template <class coeff_t>
  inline Vector<coeff_t> operator*=
  ( Vector<coeff_t>& X, const coeff_t& alpha)
  {
    X = alpha*X;
    return X;
  }

  template <class coeff_t>
  inline Vector<coeff_t> operator/=
  ( Vector<coeff_t>& X, const coeff_t& alpha)
  {
    X = X / alpha;
    return X;
  }

  template <class coeff_t>
  inline Vector<coeff_t> operator*
  (const Vector<coeff_t>& X, const coeff_t& alpha)
  { return operator*(alpha, X); }

  template <class coeff_t>
  inline Vector<coeff_t> operator/
  (const Vector<coeff_t>& X, const coeff_t& alpha)
  {
    assert(!close(alpha, static_cast<coeff_t>(0.)));
    coeff_t invalpha = static_cast<coeff_t>(1.) / alpha; 
    return operator*(invalpha, X); 
  }

  template <class coeff_t>
  inline Matrix<coeff_t> operator+
  (const Matrix<coeff_t>& X, const Matrix<coeff_t>& Y)
  {
    Matrix<coeff_t> res(Y);
    Add(X, res);
    return res;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> operator-
  (const Matrix<coeff_t>& X, const Matrix<coeff_t>& Y)
  {
    Matrix<coeff_t> res(X);
    Add(Y, res, static_cast<coeff_t>(-1.));
    return res;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> operator+
  (const Matrix<coeff_t>& X, const coeff_t& c)
  {
    Matrix<coeff_t> res(X);
    Map(res, [&c](coeff_t& x) { x = x + c; } );
    return res;
  }

  template <class coeff_t>
  inline Matrix<coeff_t>& operator +=
  (Matrix<coeff_t>& a, const Matrix<coeff_t>& b)
  {
    a = a + b;
    return a;
  }
  template <class coeff_t>
  inline Matrix<coeff_t>& operator -=
  (Matrix<coeff_t>& a, const Matrix<coeff_t>& b)
  {
    a = a - b;
    return a;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> operator-
  (const Matrix<coeff_t>& X, const coeff_t& c)
  {
    Matrix<coeff_t> res(X);
    Map(res, [&c](coeff_t& x) { x = x - c; } );
    return res;
  }


  template <class coeff_t>
  inline Matrix<coeff_t> operator-(const Matrix<coeff_t>& X)
  {
    Matrix<coeff_t> res(X);
    Scale(static_cast<coeff_t>(-1.), res);
    return res;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> operator*
  (const coeff_t& alpha,  const Matrix<coeff_t>& X)
  {
    Matrix<coeff_t> res(X);
    Scale(alpha, res);
    return res;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> operator*=
  ( Matrix<coeff_t>& X, const coeff_t& alpha)
  {
    X = alpha*X;
    return X;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> operator/=
  ( Matrix<coeff_t>& X, const coeff_t& alpha)
  {
    X = X / alpha;
    return X;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> operator*
  (const Matrix<coeff_t>& X, const coeff_t& alpha)
  { return operator*(alpha, X); }

  template <class coeff_t>
  inline Matrix<coeff_t> operator/
  (const Matrix<coeff_t>& X, const coeff_t& alpha)
  {
    assert(!close(alpha, static_cast<coeff_t>(0.)));
    coeff_t invalpha = static_cast<coeff_t>(1.) / alpha; 
    return operator*(invalpha, X); 
  }

}

#endif
