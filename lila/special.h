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
#include "range.h"

namespace lila {

  template <class matrix_t>
  inline void Zeros(matrix_t& mat)
  { std::fill(mat.begin(), mat.end(), 0.); }

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

}

#endif
