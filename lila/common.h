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

#ifndef LILA_COMMON_H_
#define LILA_COMMON_H_

#include <complex>

namespace lila {
  using int32 = int;
  using uint32 = unsigned int;
  using int64 = long;
  using uint64 = unsigned long;

  using size_type = long;

  using scomplex = std::complex<float>;
  using complex = std::complex<double>;
}

#endif
