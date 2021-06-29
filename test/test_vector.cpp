#include "catch.hpp"

#include <lila/all.h>

TEST_CASE( "vector", "[core]" ) {
  int m=5;

  lila::Vector<float> fvec(m);
  fvec(0) = 1.23;
  fvec(1) = 3.21;
  fvec(2) = 9.87;
  fvec(3) = 6.54;
  fvec(4) = 4.56;
  // Print(fvec);  
  lila::Vector<float> fvec2(fvec);  
  // Print(fvec2);
  REQUIRE(equal(fvec, fvec2));
  
  lila::Vector<double> dvec(m);
  dvec(0) = 1.23;
  dvec(1) = 3.21;
  dvec(2) = 9.87;
  dvec(3) = 6.54;
  dvec(4) = 4.56;
  // Print(dvec);  
  lila::Vector<double> dvec2(dvec);  
  REQUIRE(equal(dvec, dvec2));

  lila::Vector<std::complex<float>> cvec(m);
  cvec(0) = {1.23, 3.12};
  cvec(1) = {3.21, -1.32};
  cvec(2) = {9.87, -7.98};
  cvec(3) = {6.54, 4.65};
  cvec(4) = {4.56, 6.45};
  // Print(cvec);  
  lila::Vector<std::complex<float>> cvec2(cvec);  
  // Print(cvec2);  
  REQUIRE(equal(cvec, cvec2));

  lila::Vector<std::complex<double>> zvec(m);
  zvec(0) = {1.23, 3.12};
  zvec(1) = {3.21, -1.32};
  zvec(2) = {9.87, -7.98};
  zvec(3) = {6.54, 4.65};
  zvec(4) = {4.56, 6.45};
  // Print(zvec);  
  lila::Vector<std::complex<double>> zvec2(zvec);  
  REQUIRE(equal(zvec, zvec2));
  // Print(zvec2); 
  
  auto zvec3 = lila::String2Vector<std::complex<double>>(lila::Vector2String(zvec));
  REQUIRE(equal(zvec, zvec3));

}
