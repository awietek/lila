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

#ifndef LILA_STRINGS_H_
#define LILA_STRINGS_H_

#include <cstdlib>
#include <string>
#include <sstream>
#include <complex>

namespace lila {

  /*! @brief Converts a string tu a real/complex number.

    @param str string to be converted
    @return real/complex number generated form string
    @tparam coeff_t type of returned number 
   */
  template <class coeff_t>
  inline coeff_t string2number(const std::string& str);

  template <> inline float string2number<float>(const std::string& str)
  { return strtof(str.c_str(), NULL); }

  template <> inline double string2number<double>(const std::string& str)
  { return strtod(str.c_str(), NULL); }

  template <> inline std::complex<float> string2number<std::complex<float>>
  (const std::string& str)
  { 
    std::istringstream is(str);
    std::complex<float> c;
    is >> c;
    return c; 
  }

  template <> inline std::complex<double> string2number<std::complex<double>>
  (const std::string& str)
  { 
    std::istringstream is(str);
    std::complex<double> c;
    is >> c;
    return c; 
  }

}

#endif
