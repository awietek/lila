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

#ifndef LILA_ADD_H_
#define LILA_ADD_H_

#include <cassert>
#include <algorithm>

#include "complex.h"
#include "matrix.h"
#include "vector.h"
#include "blaslapack.h"

namespace lila {
  
  // template <class object_t, class function_t>
  // inline void Map(object_t& X, function_t func)
  // { std::for_each(X.begin(), X.end(), func); }

  template <class object_t, class function_t>
  inline object_t Map(object_t&& X, function_t func)
  { 
    std::for_each(X.begin(), X.end(), func);
    return X;
  }

  template <class coeff_t, template<class> class object_t>
  inline void Copy(const object_t<coeff_t>& X,  object_t<coeff_t>& Y)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::copy(&dx, X.data(), &inc, Y.data(), &inc);
  }

  template <class coeff_t, template<class> class object_t>
  inline void Add(const object_t<coeff_t>& X,  object_t<coeff_t>& Y, 
		  coeff_t alpha = static_cast<coeff_t>(1.))
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::axpy(&dx, &alpha, X.data(), &inc, Y.data(), &inc);
  }

  template <class coeff_t, template<class> class object_t>
  inline void Scale(const coeff_t& alpha,  object_t<coeff_t>& X)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type inc = 1;
    blaslapack::scal(&dx, &alpha, X.data(), &inc);
  }

  template <class coeff_t, template<class> class object_t>
  inline coeff_t Dot(const object_t<coeff_t>& X, const object_t<coeff_t>& Y)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    return blaslapack::dot(&dx, X.data(), &inc, Y.data(), &inc);
  }

  template <class coeff_t, template<class> class object_t>
  inline real_t<coeff_t> Norm(const object_t<coeff_t>& X)
  {
    // TODO: use proper LAPACK function here
    return sqrt(real(Dot(X,X)));
  }

  template <class coeff_t, template<class> class object_t>
  inline coeff_t Sum(const object_t<coeff_t>& X)
  { return std::accumulate(X.begin(), X.end(), 0.); }


  template <class coeff_t, template<class> class object_t>
  inline object_t<coeff_t> operator+
  (const object_t<coeff_t>& X, const object_t<coeff_t>& Y)
  {
    object_t<coeff_t> res(Y);
    Add(X, res);
    return res;
  }

  template <class coeff_t, template<class> class object_t>
  inline object_t<coeff_t> operator-
  (const object_t<coeff_t>& X, const object_t<coeff_t>& Y)
  {
    object_t<coeff_t> res(X);
    Add(Y, res, static_cast<coeff_t>(-1.));
    return res;
  }

  template <class coeff_t, template<class> class object_t>
  inline object_t<coeff_t> operator+
  (const object_t<coeff_t>& X, const coeff_t& c)
  {
    object_t<coeff_t> res(X);
    Map(res, [&c](coeff_t& x) { x = x + c; } );
    return res;
  }

  template <class coeff_t, template<class> class object_t>
  inline object_t<coeff_t> operator-
  (const object_t<coeff_t>& X, const coeff_t& c)
  {
    object_t<coeff_t> res(X);
    Map(res, [&c](coeff_t& x) { x = x - c; } );
    return res;
  }


  template <class coeff_t, template<class> class object_t>
  inline object_t<coeff_t> operator-(const object_t<coeff_t>& X)
  {
    object_t<coeff_t> res(X);
    Scale(static_cast<coeff_t>(-1.), res);
    return res;
  }

  template <class coeff_t, template<class> class object_t>
  inline object_t<coeff_t> operator*
  (const coeff_t& alpha,  const object_t<coeff_t>& X)
  {
    object_t<coeff_t> res(X);
    Scale(alpha, res);
    return res;
  }

  template <class coeff_t, template<class> class object_t>
  inline object_t<coeff_t> operator*
  (const object_t<coeff_t>& X, const coeff_t& alpha)
  { return operator*(alpha, X); }

  template <class coeff_t, template<class> class object_t>
  inline object_t<coeff_t> operator/
  (const object_t<coeff_t>& X, const coeff_t& alpha)
  {
    assert(!close(alpha, static_cast<coeff_t>(0.)));
    coeff_t invalpha = static_cast<coeff_t>(1.) / alpha; 
    return operator*(invalpha, X); 
  }
}

#endif
