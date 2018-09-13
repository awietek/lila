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

#include "complex.h"
#include "precision.h"

namespace lila {
  
  template <class coeff_t, template<class> class object_t>
  inline bool close(const object_t<coeff_t>& mat1, const object_t<coeff_t>& mat2)
  {
    assert(mat1.size() == mat2.size());
    return std::equal(mat1.begin(), mat1.end(), mat2.begin(), 
		      [](coeff_t x, coeff_t y) {
			// Floating point comparision to close zero
			return std::abs(x - y) <= 
			  precision_c<coeff_t>::val() * std::abs(x); 
		      });
  }
  
  template <class coeff_t, template<class> class object_t>
  inline bool equal(const object_t<coeff_t>& mat1, const object_t<coeff_t>& mat2)
  {
    assert(mat1.size() == mat2.size());
    return std::equal(mat1.begin(), mat1.end(), mat2.begin());
  }
  
}
#endif
