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
void test_expm()
{
  using namespace lila;
  int n=10;
  coeff_t alpha=1.23;
  for (int seed : range<int>(10)) {
  
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    Matrix<coeff_t> A(n, n);
    Random(A, fgen);
    auto Aexp = ExpM(A, alpha);

    // Check whether hermitian gives same result
    auto A1 = A;
    Add(Herm(A1), A1);

    // auto LILA_CLK(t1);
    Aexp = ExpM(A1, alpha);
    // auto LILA_CLK(t2);
    // LILA_TIME_MICRO("pade", t1, t2);

    auto Aexp2 = A1;
    // LILA_CLK(t1);
    ExpSym(Aexp2, alpha);
    // LILA_CLK(t2);
    // LILA_TIME_MICRO("herm", t1, t2);
    // LilaPrint(Aexp);
    // LilaPrint(Aexp2);
    REQUIRE(close(Aexp2, Aexp));

    // Pade from EXPOKIT is actually really good in timing!!
    
  }
}


TEST_CASE( "Matrix ExpM test", "[expm]" ) {
  // test_expm<float>();
  test_expm<double>();
  // test_expm<std::complex<float>>();
  test_expm<std::complex<double>>();
}
