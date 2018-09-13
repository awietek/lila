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

#ifndef LILA_PRECISION_H_
#define LILA_PRECISION_H_

#include "complex.h"

namespace lila {
  
  // Predifined relative precision for comparisions
  static constexpr float float_precision_c = 1e-4;
  static constexpr double double_precision_c = 1e-12;

  template<class T> 
  struct precision_c {
    static constexpr lila::real_t<T> val();
  };

  template<> 
  struct precision_c<float> {
    static constexpr lila::real_t<float> val() 
    { return float_precision_c; }
  };

  template<> 
  struct precision_c<double> {
    static constexpr lila::real_t<double> val() 
    { return double_precision_c; }
  };
  
  template<> 
  struct precision_c<std::complex<float>> {
    static constexpr lila::real_t<std::complex<float>> val() 
    { return float_precision_c; }
  };

  template<> 
  struct precision_c<std::complex<double>> {
    static constexpr lila::real_t<std::complex<double>> val() 
    { return double_precision_c; }
  };
  

}

#endif
