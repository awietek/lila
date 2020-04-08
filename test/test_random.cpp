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
#include <lila/print.h>
#include <lila/mult.h>
#include <lila/compare.h>

TEST_CASE( "Random Matrix test", "[Random]" ) {

  int m=3;
  int n=4;

  int seed=42;
  
  lila::uniform_dist_t<float> fdist(-1., 1.);
  lila::uniform_gen_t<float> fgen(fdist, seed);
  lila::uniform_gen_t<std::complex<float>> cgen(fdist, seed);
  
  lila::uniform_dist_t<double> ddist(-1., 1.);
  lila::uniform_gen_t<double> dgen(ddist, seed);
  lila::uniform_gen_t<std::complex<double>> zgen(ddist, seed);
  
  lila::Matrix<float> fmat1(m, n);
  Random(fmat1, fgen);
  // Print(fmat1);

  lila::Matrix<double> dmat1(m, n);
  Random(dmat1, dgen);
  // Print(dmat1);

  lila::Matrix<std::complex<float>> cmat1(m, n);
  Random(cmat1, cgen);
  // Print(cmat1);

  lila::Matrix<std::complex<double>> zmat1(m, n);
  Random(zmat1, zgen);
  // Print(zmat1);
  
  // reset seed
  fgen.seed(seed);
  dgen.seed(seed);

  lila::Matrix<float> fmat2(m, n);
  Random(fmat2, fgen);
  // Print(fmat2);

  lila::Matrix<double> dmat2(m, n);
  Random(dmat2, dgen);
  // Print(dmat2);

  lila::Matrix<std::complex<float>> cmat2(m, n);
  Random(cmat2, cgen);
  // Print(cmat2);

  lila::Matrix<std::complex<double>> zmat2(m, n);
  Random(zmat2, zgen);
  // Print(zmat2);

  m=5; 
  n=5;
  // {
  //   lila::Matrix<float> umat(m, n);
  //   lila::RandomUnitary(umat, fgen);
  //   lila::Matrix<float> p(m, m);
  //   lila::Mult(umat, umat, p, float(1.), float(0.), 'C', 'N');
  //   lila::Matrix<float> Id(p);
  //   lila::Identity(Id);
  //   REQUIRE(lila::close<float>(Id, p));
  // }
  {
    lila::Matrix<double> umat(m, n);
    lila::RandomUnitary(umat, fgen);
    lila::Matrix<double> p(m, m);
    lila::Mult(umat, umat, p, 1., 0., 'C', 'N');
    lila::Matrix<double> Id(p);
    lila::Identity(Id);
    REQUIRE(lila::close<double>(Id, p));
  }
  // {
  //   lila::Matrix<std::complex<float>> umat(m, n);
  //   lila::RandomUnitary(umat, cgen);
  //   lila::Matrix<std::complex<float>> p(m, m);
  //   lila::Mult(umat, umat, p, std::complex<float>(1.), 
  // 	       std::complex<float>(0.), 'C', 'N');
  //   lila::Matrix<std::complex<float>> Id(p);
  //   lila::Identity(Id);
  //   REQUIRE(lila::close<std::complex<float>>(Id, p));
  // }
  {
    lila::Matrix<std::complex<double>> umat(m, n);
    lila::RandomUnitary(umat, zgen);
    lila::Matrix<std::complex<double>> p(m, m);
    lila::Mult(umat, umat, p, std::complex<double>(1.), 
	       std::complex<double>(0.), 'C', 'N');
    lila::Matrix<std::complex<double>> Id(p);
    lila::Identity(Id);
    REQUIRE(lila::close<std::complex<double>>(Id, p));
  }
}
