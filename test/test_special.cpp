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
#include <lila/special.h>
#include <lila/print.h>
#include <lila/range.h>

template <class matrix_t>
bool all_zero(const matrix_t& mat)
{
  bool all_zero=true;
  for (auto i : mat.rows())
    for (auto j : mat.cols())
      all_zero &= (mat(i,j) == static_cast<typename matrix_t::coeff_type>(0.));
  return all_zero;
}

template <class matrix_t>
bool all_one(const matrix_t& mat)
{
  bool all_zero=true;
  for (auto i : mat.rows())
    for (auto j : mat.cols())
      all_zero &= (mat(i,j) == static_cast<typename matrix_t::coeff_type>(1.));
  return all_zero;
}

TEST_CASE( "Special Matrix test", "[Special]" ) {

  int m=3;
  int n=4;
  
  // Zeros
  {
    lila::Matrix<float> fmat(m, n);
    Zeros(fmat);
    REQUIRE(all_zero(fmat));

    lila::Matrix<double> dmat(m, n);
    Zeros(dmat);
    REQUIRE(all_zero(dmat));

    lila::Matrix<std::complex<float>> cmat(m, n);
    Zeros(cmat);
    REQUIRE(all_zero(cmat));

    lila::Matrix<std::complex<double>> zmat(m, n);
    Zeros(zmat);
    REQUIRE(all_zero(zmat));
  }

  // Ones
  {
    lila::Matrix<float> fmat(m, n);
    Ones(fmat);
    REQUIRE(all_one(fmat));

    lila::Matrix<double> dmat(m, n);
    Ones(dmat);
    REQUIRE(all_one(dmat));

    lila::Matrix<std::complex<float>> cmat(m, n);
    Ones(cmat);
    REQUIRE(all_one(cmat));

    lila::Matrix<std::complex<double>> zmat(m, n);
    Ones(zmat);
    REQUIRE(all_one(zmat));
  }
  
  // Unitary
  int N = 3;
  lila::Vector<double> params(N*N);
  for (int i : lila::range<int>(N*N))
    params(i) = 0.1 * (i+1);
  auto unitary = lila::Unitary<std::complex<double>>(N, params);
  // LilaPrint(unitary);
  REQUIRE(IsUnitary(unitary));

}
