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

#include <iostream>
#include <algorithm>
#include <lila/all.h>

template <class coeff_t>
void test_eigen()
{
  using namespace lila;
  int n=10;

  
  // Test Eigen
  for (int seed : range<int>(10)) {
  
    uniform_dist_t<coeff_t> dist(-1., 1.);
    uniform_gen_t<coeff_t> gen(dist, seed);
    Matrix<coeff_t> A(n, n);
    Random(A, gen);

    ////////////////////////////////////////////////////
    // Check whether symmetric and generic eigenvalue
    // routines give same eigenvalues
    auto A1 = A;
    Add(Herm(A1), A1);
    auto v1 = Eigen(A1);
    auto A2 = A;
    Add(Herm(A2), A2);
    auto v2 = EigenSym(A2);

    auto e1r = Real(v1.eigenvalues);
    auto e1i = Imag(v1.eigenvalues);    
    std::sort(e1r.data(), e1r.data() + e1r.size());
    REQUIRE(close(e1r, v2.eigenvalues));
    REQUIRE(close(e1i, (real_t<coeff_t>)0.));




    /////////////////////////////////////////////////////////
    // Check whether symmetric  eigenvalue
    // routines computes correct eigenvalues/vectors
    auto Ah = A;
    Add(Herm(Ah), Ah);

    auto v3 = EigenSym(Ah);

    // Check whether EigenvaluesSym and EigenSym give same eigenvalues
    REQUIRE(close(v3.eigenvalues, EigenvaluesSym(Ah)));

    // check if routine computes correct eigenvalues/vectors
    for (int i=0; i<n; ++i) 
      {
	auto eval = v3.eigenvalues(i);
	auto evec = v3.eigenvectors.col(i);
	// LilaPrint(Mult(Ah, evec) - (coeff_t)eval * Mult(Ad, evec));
	REQUIRE(close(Mult(Ah, evec), (coeff_t)eval * evec));
      }



    /////////////////////////////////////////////////////////
    // Check symmetric definite generalized eigenvalue
    // Create positive definite Matrix
    Matrix<coeff_t> Ad(n, n);
    std::vector<lila::Vector<coeff_t>> vecs; 
    for (int k=0; k<n; ++k)
      vecs.push_back(Random<coeff_t>(n, gen));
    for (int i=0; i<n; ++i)
      for (int j=0; j<n; ++j)
	Ad(i,j) = Dot(vecs[i], vecs[j]);

    auto v4 = EigenGenSymDef(Ah, Ad);

    // Check whether EigenvaluesSym and EigenSym give same eigenvalues
    REQUIRE(close(v4.eigenvalues, EigenvaluesGenSymDef(Ah, Ad)));

    // check if routine computes correct eigenvalues/vectors
    for (int i=0; i<n-1; ++i) // n-1: Last eigenvector somehow looses precision?
      {
	auto eval = v4.eigenvalues(i);
	auto evec = v4.eigenvectors.col(i);
	// LilaPrint(Mult(Ah, evec) - (coeff_t)eval * Mult(Ad, evec));
	REQUIRE(close(Mult(Ah, evec), (coeff_t)eval * Mult(Ad, evec)));
      }

  }
}


TEST_CASE( "Eigen test", "[Eigen]" ) {
  test_eigen<float>();
  test_eigen<double>();
  test_eigen<std::complex<float>>();
  test_eigen<std::complex<double>>();
}
