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

#ifndef LILA_COMPLEX_H_
#define LILA_COMPLEX_H_

#include <complex>

#include "detail/complex_detail.h"

namespace lila {
  
  template <class coeff_t>
  using real_t = typename detail::real_type_struct<coeff_t>::type;

  template <class T> T real (T x) { return x; }
  template <class T> T real (std::complex<T> x) { return x.real(); }
  template <class T> T conj (T x) { return x; }
  template <class T> std::complex<T> conj(std::complex<T> x)
  { return std::conj(x); }

}
#endif  // LILA_COMPLEX_H_
