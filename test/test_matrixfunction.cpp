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

#include <lila/matrixfunction.h>
#include <lila/print.h>
#include <lila/random.h>
#include <lila/add.h>
#include <lila/compare.h>

template <class coeff_t>
void test_matrixfunction()
{
  using namespace lila;
  int n=5;
  coeff_t alpha=1.23;
  for (int seed : range<int>(10)) {
  
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    Matrix<coeff_t> A(n, n);
    Random(A, fgen);
    auto A1 = A;
    Add(Herm(A1), A1);
    auto Aexp = A1;
    ExpSym(Aexp, alpha);
    auto Aret = Aexp;
    LogSym(Aret, alpha);

    // LilaPrint(A1);
    // LilaPrint(Aret);
    REQUIRE(close(A1, Aret));

    auto v2 = EigenvaluesSym(Aexp);
    for (auto e : v2)
      REQUIRE(e >= 0);
    
  }
}


TEST_CASE( "Matrix function test", "[matrixfunction]" ) {
  // test_matrixfunction<float>();
  test_matrixfunction<double>();
  // test_matrixfunction<std::complex<float>>();
  test_matrixfunction<std::complex<double>>();
}
