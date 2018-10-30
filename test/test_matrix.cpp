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
#include <lila/compare.h>
#include <lila/print.h>

TEST_CASE( "Basic Matrix test", "[Matrix]" ) {
  int m=2;
  int n=3;

  lila::Matrix<float> fmat(m, n);
  fmat(0,0) = 1.23;
  fmat(0,1) = 3.21;
  fmat(0,2) = 9.87;
  fmat(1,0) = 6.54;
  fmat(1,1) = 4.56;
  fmat(1,2) = 7.89;     
  // LilaPrint(fmat);  
  // LilaPrint(Transpose(fmat));  
  // LilaPrint(Conj(fmat));  
  // LilaPrint(Herm(fmat));  
  lila::Matrix<float> fmat2(fmat);  
  REQUIRE(close(fmat, fmat2));
  
  lila::Matrix<double> dmat(m, n);
  dmat(0,0) = 1.23;
  dmat(0,1) = 3.21;
  dmat(0,2) = 9.87;
  dmat(1,0) = 6.54;
  dmat(1,1) = 4.56;
  dmat(1,2) = 7.89;     
  // LilaPrint(dmat);  
  // LilaPrint(Transpose(dmat));  
  // LilaPrint(Conj(dmat));  
  // LilaPrint(Herm(dmat));  
  lila::Matrix<double> dmat2(dmat);  
  REQUIRE(close(dmat, dmat2));

  lila::Matrix<std::complex<float>> cmat(m, n);
  cmat(0,0) = {1.23, 3.12};
  cmat(0,1) = {3.21, -1.32};
  cmat(0,2) = {9.87, -7.98};
  cmat(1,0) = {6.54, 4.65};
  cmat(1,1) = {4.56, 6.45};
  cmat(1,2) = {7.89, 9.78};     
  // LilaPrint(cmat);  
  // LilaPrint(Transpose(cmat));  
  // LilaPrint(Conj(cmat));  
  // LilaPrint(Herm(cmat));   
  lila::Matrix<std::complex<float>> cmat2(cmat);  
  REQUIRE(close(cmat, cmat2));

  lila::Matrix<std::complex<double>> zmat(m, n);
  zmat(0,0) = {1.23, 3.12};
  zmat(0,1) = {3.21, -1.32};
  zmat(0,2) = {9.87, -7.98};
  zmat(1,0) = {6.54, 4.65};
  zmat(1,1) = {4.56, 6.45};
  zmat(1,2) = {7.89, 9.78};     
  // LilaPrint(zmat);  
  // LilaPrint(Transpose(zmat));  
  // LilaPrint(Conj(zmat));  
  // LilaPrint(Herm(zmat));  
  lila::Matrix<std::complex<double>> zmat2(zmat);  
  REQUIRE(close(zmat, zmat2));
}
