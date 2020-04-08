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

#include <iostream>
template <class coeff_t>
void test_cholesky()
{
  using namespace lila;
  int n=10;
  
  // Test Cholesky
  for (int seed : range<int>(10)) 
    {
      Matrix<coeff_t> A(n, n);
    
      // Create positive definite Matrix
      uniform_dist_t<coeff_t> dist(-1., 1.);
      uniform_gen_t<coeff_t> gen(dist, seed);
      std::vector<lila::Vector<coeff_t>> vecs; 
      for (int k=0; k<n; ++k)
	vecs.push_back(Random<coeff_t>(n, gen));
      for (int i=0; i<n; ++i)
	for (int j=0; j<n; ++j)
	  A(i,j) = Dot(vecs[i], vecs[j]);
    
      // LilaPrint(A);
      auto U = Cholesky(A, 'U');
      // LilaPrint(U);
      auto produ = Mult(Herm(U), U);
      // LilaPrint(produ);
      REQUIRE(close(produ, A));      

      auto L = Cholesky(A, 'L');
      // LilaPrint(L);
      auto prodl = Mult(L, Herm(L));
      // LilaPrint(prodl);
      REQUIRE(close(prodl, A)); 

    }
}


TEST_CASE( "Cholesky test", "[Cholesky]" ) {
  // test_cholesky<float>();
  test_cholesky<double>();
  // test_cholesky<std::complex<float>>();
  test_cholesky<std::complex<double>>();
}
