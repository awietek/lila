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
void test_lanczos()
{
  using namespace lila;

  int n=100;
  double precision = 1e-4;
  int max_iterations = n;
  int n_lowest = 3;
  int num_eigenvalue = n_lowest;

  // Test Lanczos
  for (int random_seed : range<int>(10)) {
  
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, random_seed);
  
    Matrix<coeff_t> A(n, n);
    Random(A, fgen);
    A += Herm(A);

    REQUIRE(close(A, Herm(A)));
    auto true_eigs = EigenvaluesH(A);

    
    uint64 dim = A.nrows();
    auto multiply = [&A](const Vector<coeff_t>& v, Vector<coeff_t>& w) {
      w = Mult(A, v);
    };
    auto lzs = Lanczos<coeff_t, decltype(multiply)>
      (dim, random_seed, max_iterations, precision, num_eigenvalue, multiply);
    auto lzs_eigenvalues = lzs.eigenvalues();
   
    std::vector<int> num_eigenvectors;
    for (int k = 0; k<n_lowest; ++k)
      num_eigenvectors.push_back(k);

    // LilaPrint(lzs.alphas().size());
    auto eigenvectors = lzs.eigenvectors(num_eigenvectors);
    for (int k = 0; k<n_lowest; ++k)
      {
	// LilaPrint(k);
	// LilaPrint(random_seed);
	// LilaPrint(true_eigs(k));
	// LilaPrint(lzs_eigenvalues(k));
	// LilaPrint(Dot(eigenvectors[k], Mult(A, eigenvectors[k])));
	// LilaPrint(Norm(Mult(A, eigenvectors[k]) - (coeff_t)lzs_eigenvalues(k)*eigenvectors[k]));
	// REQUIRE(Norm(Mult(A, eigenvectors[k]) - lzs_eigenvalues(k)*eigenvectors[k]) < precision);
	REQUIRE(std::abs(true_eigs(k) - lzs_eigenvalues(k)) < 10*precision);
      }
  }
}


TEST_CASE( "Lanczos test", "[Lanczos]" ) {
  // test_lanczos<float>();
  test_lanczos<double>();
  // test_lanczos<std::complex<float>>();
  test_lanczos<std::complex<double>>();
}
