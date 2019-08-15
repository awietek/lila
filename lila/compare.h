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

#ifndef LILA_COMPARE_H_
#define LILA_COMPARE_H_

#include <cmath>  
#include <cassert>

#include "matrix.h"
#include "vector.h"
#include "complex.h"
#include "precision.h"

namespace lila {
  
  template <class coeff_t>
  inline bool close(const coeff_t& x, const coeff_t& y,
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  { return (std::abs(x - y) <= ( atol + rtol * std::abs(y))); }
  
  template <class coeff_t>
  inline bool close(const Matrix<coeff_t>& mat1, const Matrix<coeff_t>& mat2, 
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  {
    assert(mat1.nrows() == mat2.nrows());
    assert(mat1.ncols() == mat2.ncols());
    return std::equal(mat1.begin(), mat1.end(), mat2.begin(),
		      [&atol, &rtol](coeff_t x, coeff_t y) {
		      	return close(x, y, atol, rtol); 
		      });
  }

  template <class coeff_t>
  inline bool close(const Vector<coeff_t>& vec1, const Vector<coeff_t>& vec2, 
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  {
    assert(vec1.size() == vec2.size());
    return std::equal(vec1.begin(), vec1.end(), vec2.begin(),
		      [&atol, &rtol](coeff_t x, coeff_t y) {
		      	return close(x, y, atol, rtol); 
		      });
  }
  
  template <class coeff_t>
  inline bool close(const Matrix<coeff_t>& mat, const coeff_t& val, 
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  {
    return std::all_of(mat.begin(), mat.end(), 
		       [&atol, &rtol, &val](coeff_t x) {
		      	return close(x, val, atol, rtol); 
		      });
  }

  template <class coeff_t>
  inline bool close(const Vector<coeff_t>& vec, const coeff_t& val, 
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  {
    return std::all_of(vec.begin(), vec.end(),
		       [&atol, &rtol, &val](coeff_t x) {
			 return close(x, val, atol, rtol); 
		       });
  }


  template <class coeff_t, template<class> class object_t>
  inline bool equal(const object_t<coeff_t>& mat1,
		    const object_t<coeff_t>& mat2)
  {
    assert(mat1.size() == mat2.size());
    return std::equal(mat1.begin(), mat1.end(), mat2.begin());
  }
  
}
#endif
