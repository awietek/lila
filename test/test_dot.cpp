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


#include "catch.hpp"

#include <lila/all.h>

template <class coeff_t>
void test_dot()
{
  using namespace lila;
  int n=20;
  
  // Test Dot
  for (int seed : lila::range<int>(10)) {
  
    lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
    lila::uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    lila::Vector<coeff_t> v1(n);
    Random(v1, fgen);
   
    lila::Vector<coeff_t> v2(n);
    Random(v2, fgen);

    coeff_t dot_lila = Dot(v1, v2);
    coeff_t dot_test = 0;
    for (int i = 0; i < n; ++i)
      dot_test += lila::conj(v1(i)) * v2(i);

    REQUIRE(close(dot_lila, dot_test));
  }
}


TEST_CASE( "Dot test", "[Dot]" ) {
  // test_dot<float>();
  test_dot<double>();
  // test_dot<std::complex<float>>();
  test_dot<std::complex<double>>();
}
