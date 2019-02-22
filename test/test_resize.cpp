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

template <class coeff_t>
void test_resize()
{
  using namespace lila;

  int m=6;
  int n=5;

  // Test Resize
  for (int m_new=0; m_new < 10; ++m_new)
    for (int n_new=0; n_new < 10; ++n_new)
      {
	for (int seed : lila::range<int>(3)) 
	  {
      
	    uniform_dist_t<coeff_t> dist(-1., 1.);
	    uniform_gen_t<coeff_t> gen(dist, seed);
      
	    Matrix<coeff_t> A(m, n);
	    Random(A, gen);
	    auto A_res = A;
	    A_res.resize(m_new, n_new);
	    for (int i=0; i < m_new; ++i)
	      for (int j=0; j < n_new; ++j)
		{
		  if ((i < std::min(m_new, m)) && (j < std::min(n_new, n))) 
		    REQUIRE(A_res(i, j) == A(i, j));
		  else REQUIRE(close<coeff_t>(A_res(i, j), 0.));
		}
      	  }
      }
}


TEST_CASE( "Resize test", "[Resize]" ) {
  test_resize<float>();
  test_resize<double>();
  test_resize<std::complex<float>>();
  test_resize<std::complex<double>>();
}
