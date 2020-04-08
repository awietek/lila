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

#include <string>

#include <lila/matrix.h>
#include <lila/random.h>
#include <lila/special.h>
#include <lila/mult.h>
#include <lila/print.h>
#include <lila/compare.h>


template <class coeff_t>
void test_string()
{
  int m=7;
  int k=9;

  int seed=42;
  
  lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
  lila::uniform_gen_t<coeff_t> fgen(fdist, seed);
  lila::Matrix<coeff_t> A(m, k);
  Random(A, fgen);
  // Print(A);
  std::string str = lila::WriteMatrix(A);
  lila::Matrix<coeff_t> B = lila::ParseMatrix<coeff_t>(str);
  // Print(B);
}


TEST_CASE( "Matrix string conversion test", "[String]" ) {
  // test_string<float>();
  test_string<double>();
  // test_string<std::complex<float>>();
  test_string<std::complex<double>>();
}
