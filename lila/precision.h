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

#include <limits>
#include "complex.h"

namespace lila {
  
  // Predifined absolute precision for comparisions
  static constexpr float float_atol = 
    1000.*std::numeric_limits<float>::epsilon();
  static constexpr double double_atol = 
    1000.*std::numeric_limits<double>::epsilon();

  template<class T> 
  struct atol {
    static constexpr lila::real_t<T> val();
  };

  template<> 
  struct atol<float> {
    static constexpr lila::real_t<float> val() 
    { return float_atol; }
  };

  template<> 
  struct atol<double> {
    static constexpr lila::real_t<double> val() 
    { return double_atol; }
  };
  
  template<> 
  struct atol<std::complex<float>> {
    static constexpr lila::real_t<std::complex<float>> val() 
    { return float_atol; }
  };

  template<> 
  struct atol<std::complex<double>> {
    static constexpr lila::real_t<std::complex<double>> val() 
    { return double_atol; }
  };
 
  // Predifined relative precision for comparisions
  static constexpr float float_rtol = 
    10000.*std::numeric_limits<float>::epsilon();
  static constexpr double double_rtol = 
    10000.*std::numeric_limits<double>::epsilon();

  template<class T> 
  struct rtol {
    static constexpr lila::real_t<T> val();
  };

  template<> 
  struct rtol<float> {
    static constexpr lila::real_t<float> val() 
    { return float_rtol; }
  };

  template<> 
  struct rtol<double> {
    static constexpr lila::real_t<double> val() 
    { return double_rtol; }
  };
  
  template<> 
  struct rtol<std::complex<float>> {
    static constexpr lila::real_t<std::complex<float>> val() 
    { return float_rtol; }
  };

  template<> 
  struct rtol<std::complex<double>> {
    static constexpr lila::real_t<std::complex<double>> val() 
    { return double_rtol; }
  };
  
  

}

#endif
