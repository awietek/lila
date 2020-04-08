// Copyright 2019 Alexander Wietek - All Rights Reserved.
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
#include <iostream>

template <class coeff_t>
void test_tmatrix()
{
  using namespace lila;

  int n=10;
  
  // Test Eigen
  for (int seed : lila::range<int>(10)) {
  
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    Vector<coeff_t> alphas(n);
    Vector<coeff_t> betas(n-1);
    Random(alphas, fgen);
    Random(betas, fgen);

    Tmatrix<coeff_t> tmat(alphas, betas);
    auto fullmat = Matrix<coeff_t>(tmat);
    auto eigs = Eigenvalues(tmat);
    auto res = Eigen(tmat);

    // LilaPrint(tmat.alphas());
    // LilaPrint(tmat.betas());
    // LilaPrint(fullmat);
    // LilaPrint(eigs);

    auto evecs = res.eigenvectors;
    for (int i : evecs.cols())
      {
	// printf("norm: %f\n", Dot(evecs.col(i), evecs.col(i)));
	// printf("eig : %f\n", eigs(i));
	coeff_t test = Dot(evecs.col(i), Mult(fullmat, evecs.col(i)));
	// printf("test: %f\n", test);
	REQUIRE(close(eigs(i), test));
      }

  }
}


TEST_CASE( "Tmatrix test", "[Tmatrix]" ) {
  // test_tmatrix<float>();
  test_tmatrix<double>();
}
