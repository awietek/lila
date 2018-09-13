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
#include <lila/mult.h>
#include <lila/print.h>
#include <lila/compare.h>

template <class coeff_t>
lila::Matrix<coeff_t> simple_multiply(const lila::Matrix<coeff_t>& A, 
				      const lila::Matrix<coeff_t>& B)
{
  lila::Matrix<coeff_t> C(A.nrows(), B.ncols());
  lila::Zeros(C);  
  for (auto i : C.rows())
    for (auto j : C.cols())
      for (auto k : A.cols())
	C(i,j) += A(i,k)*B(k,j);
  return C;
}

template <class coeff_t>
void test_mult()
{
  int m=7;
  int k=9;
  int n=11;

  int seed=42;
  
  lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
  lila::uniform_gen_t<coeff_t> fgen(fdist, seed);
  
  lila::Matrix<coeff_t> A(m, k);
  lila::Matrix<coeff_t> B(k, n);
  lila::Matrix<coeff_t> C(m, n);
  Random(A, fgen);  
  Random(B, fgen);
  Random(C, fgen);
  Mult(A, B, C);
  // Print(A);
  // Print(B);
  // Print(C);
  lila::Matrix<coeff_t> C2 = simple_multiply(A, B);
  // Print(C2);
  REQUIRE(lila::close<coeff_t>(C, C2));

  Mult(B, A, C, coeff_t(.5), coeff_t(.6), 'T', 'T');
  // Print(A);
  // Print(B);
  // Print(C);

  A.resize(m, k);
  B.resize(n, k);
  Mult(A, B, C, coeff_t(.5), coeff_t(.6), 'N', 'T');
  // Print(A);
  // Print(B);
  // Print(C);

  A.resize(k, m);
  B.resize(k, n);
  Mult(A, B, C, coeff_t(.5), coeff_t(.6), 'T', 'N');
  // Print(A);
  // Print(B);
  // Print(C);
}


TEST_CASE( "Matrix multiplication test", "[Mult]" ) {
  test_mult<float>();
  test_mult<double>();
  test_mult<std::complex<float>>();
  test_mult<std::complex<double>>();
}
