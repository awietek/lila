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

#include <lila/matrix.h>
#include <lila/random.h>
#include <lila/special.h>
#include <lila/eigen.h>
#include <lila/print.h>
#include <lila/compare.h>
#include <lila/range.h>
#include <lila/add.h>

#include <iostream>
template <class coeff_t>
void test_eigen()
{
  int n=5;
  
  // Test Eigen
  for (int seed : lila::range<int>(10)) {
  
    lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
    lila::uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    lila::Matrix<coeff_t> A(n, n);
    Random(A, fgen);
    auto A1 = A;
    Add(Herm(A1), A1);
    auto v1 = lila::Eigen(A1);
    auto A2 = A;
    Add(Herm(A2), A2);
    auto v2 = lila::EigenH(A2);

    auto A3 = A;
    Add(Herm(A3), A3);
    auto v3 = lila::Eigenvalues(A3);
    auto A4 = A;
    Add(Herm(A4), A4);
    auto v4 = lila::EigenvaluesH(A4);
    
    LilaPrint(Real(v1.eigenvalues));
    LilaPrint(v2);
    LilaPrint(Real(v3));
    LilaPrint(v4);
  }
}


TEST_CASE( "Eigen test", "[Eigen]" ) {
  // test_eigen<float>();
  // test_eigen<double>();
  test_eigen<std::complex<float>>();
  test_eigen<std::complex<double>>();
}
