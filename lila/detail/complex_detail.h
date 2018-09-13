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

#ifndef LILA_DETAIL_COMPLEX_DETAIL_H_
#define LILA_DETAIL_COMPLEX_DETAIL_H_

#include <complex>

namespace lila {
  namespace detail {
    
    template <class coeff_t>
    struct real_type_struct {
      typedef coeff_t type;
    };

    template <class coeff_t>
    struct real_type_struct<std::complex<coeff_t> > {
      typedef coeff_t type;
    };

  }

}
#endif  // LILA_DETAIL_COMPLEX_DETAIL_H_
