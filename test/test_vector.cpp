#include "catch.hpp"

#include <lila/all.h>

TEST_CASE( "vector", "[core]" ) {
  lila::Log("Test vector");

  int m=5;

  lila::Vector<float> fvec(m);
  fvec(0) = 1.23;
  fvec(1) = 3.21;
  fvec(2) = 9.87;
  fvec(3) = 6.54;
  fvec(4) = 4.56;
  lila::Vector<float> fvec2(fvec);  
  REQUIRE(equal(fvec, fvec2));
  auto fvec3 = fvec;
  REQUIRE(equal(fvec, fvec3));
  
  lila::Vector<double> dvec(m);
  dvec(0) = 1.23;
  dvec(1) = 3.21;
  dvec(2) = 9.87;
  dvec(3) = 6.54;
  dvec(4) = 4.56;
  lila::Vector<double> dvec2(dvec);  
  REQUIRE(equal(dvec, dvec2));
  auto dvec3 = dvec;
  REQUIRE(equal(dvec, dvec3));

  lila::Vector<std::complex<float>> cvec(m);
  cvec(0) = {1.23, 3.12};
  cvec(1) = {3.21, -1.32};
  cvec(2) = {9.87, -7.98};
  cvec(3) = {6.54, 4.65};
  cvec(4) = {4.56, 6.45};
  lila::Vector<std::complex<float>> cvec2(cvec);  
  REQUIRE(equal(cvec, cvec2));
  auto cvec3 = cvec;
  REQUIRE(equal(cvec, cvec3));

  lila::Vector<std::complex<double>> zvec(m);
  zvec(0) = {1.23, 3.12};
  zvec(1) = {3.21, -1.32};
  zvec(2) = {9.87, -7.98};
  zvec(3) = {6.54, 4.65};
  zvec(4) = {4.56, 6.45};
  lila::Vector<std::complex<double>> zvec2(zvec);  
  REQUIRE(equal(zvec, zvec2));
  auto zvec3 = zvec;
  REQUIRE(equal(zvec, zvec3));
}
