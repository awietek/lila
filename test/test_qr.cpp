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
#include <lila/solve.h>
#include <lila/print.h>
#include <lila/compare.h>
#include <lila/range.h>

template <class coeff_t>
void test_qr(int m, int n)
{
 
  // Test direct solve
  for (int seed : lila::range<int>(10)) {
  
    lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
    lila::uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    lila::Matrix<coeff_t> A(m, n);
    lila::Matrix<coeff_t> Q(m, n);
    lila::Matrix<coeff_t> R(m, n);

    Random(A, fgen);

    lila::Matrix<coeff_t> A_save(A);
    
    std::vector<coeff_t> tau = QRDecompose(A);
    Q = QRGetQ(A, tau);
    R = GetUpper(A);
    
    // Check if valid decomposition
    lila::Matrix<coeff_t> Prod(m, n);
    Mult(Q, R, Prod);
    REQUIRE(lila::close<coeff_t>(A_save, Prod));

    // check if R is upper trapezoidal
    int k = std::min(m, n);
    for (int row=0; row < k; ++row)
      for (int column=0; column < n; ++column)
	if (row > column) REQUIRE(R(row, column) == coeff_t(0.));
    
    lila::Matrix<coeff_t> Prod2(m, m);
    Mult(Q, Q, Prod2, coeff_t(1.), coeff_t(0.), 'C', 'N');

    lila::Matrix<coeff_t> Prod3(k, k);
    Mult(Q, Q, Prod3, coeff_t(1.), coeff_t(0.), 'N', 'C');

 
    lila::Matrix<coeff_t> Id2(Prod2);
    Identity(Id2);
    REQUIRE(lila::close<coeff_t>(Id2, Prod2));
    if (m <= n)
      {
	lila::Matrix<coeff_t> Id3(Prod3);
	Identity(Id3);
	REQUIRE(lila::close<coeff_t>(Id3, Prod3));
      }
  }
  
}


TEST_CASE( "QR Decomposition test", "[QR]" ) {
  // test_qr<float>(6,5);
  test_qr<double>(6,5);
  // test_qr<std::complex<float>>(6,5);
  test_qr<std::complex<double>>(6,5);
  // test_qr<float>(6,6);
  test_qr<double>(6,6);
  // test_qr<std::complex<float>>(6,6);
  test_qr<std::complex<double>>(6,6);
  // test_qr<float>(5,6);
  test_qr<double>(5,6);
  // test_qr<std::complex<float>>(5,6);
  test_qr<std::complex<double>>(5,6);
}
