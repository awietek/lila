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

#ifndef LILA_BLASLAPACK_TYPES_H_
#define LILA_BLASLAPACK_TYPES_H_
    
#include <complex>
#include "../common.h"

namespace lila { 
  namespace blaslapack {

#ifdef LILA_USE_MKL // Using the Intel MKL
#include "mkl.h"
    using blas_size_t = int;
    using blas_float_t = float;
    using blas_double_t = double;
    using blas_scomplex_t = MKL_Complex8;
    using blas_complex_t = MKL_Complex16;
#else  // Using normal LAPACK  
    using blas_size_t = int;
    using blas_float_t = float;
    using blas_double_t = double;
    using blas_scomplex_t = std::complex<float>;
    using blas_complex_t = std::complex<double>;
#endif
  }
}

#endif
