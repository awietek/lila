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

#include "blaslapack.h"

namespace lila {

  template <class coeff_t, template<class> class object_t>
  inline void copy(const object_t<coeff_t>& X,  object_t<coeff_t>& Y)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::copy(&dx, X.data(), &inc, X.data(), &inc);
  }

  template <class coeff_t, template<class> class object_t>
  inline void add(const object_t<coeff_t>& X,  object_t<coeff_t>& Y, 
		  coeff_t alpha = static_cast<coeff_t>(1.))
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::axpy(&dx, &alpha, X.data(), &inc, X.data(), &inc);
  }

  template <class coeff_t, template<class> class object_t>
  inline void scale(const coeff_t& alpha,  object_t<coeff_t>& X)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type inc = 1;
    blaslapack::scal(&dx, &alpha, X.data(), &inc);
  }

  template <class coeff_t, template<class> class object_t>
  inline void dot(const object_t<coeff_t>& X, const object_t<coeff_t>& Y)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = X.size();
    const size_type dy = y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    return blaslapack::scal(&dx, X.data(), &inc, Y.data(), &inc);
  }

  template <class coeff_t, template<class> class object_t>
  inline void norm(const object_t<coeff_t>& X)
  {
    // TODO: use proper LAPACK function here
    return sqrt(real(blaslapack::dot(X,X)));
  }

}

#endif
