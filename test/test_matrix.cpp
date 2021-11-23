#include "catch.hpp"

#include <lila/all.h>

TEST_CASE( "matrix", "[core]" ) {
  lila::Log("Test matrix");

  int m=4;
  int n=5;

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
  auto fmat3 = fmat;
  REQUIRE(close(fmat, fmat3));

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
  auto dmat3 = dmat;
  REQUIRE(close(dmat, dmat3));


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
  auto cmat3 = cmat;
  REQUIRE(close(cmat, cmat3));


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
  auto zmat3 = zmat;
  REQUIRE(close(zmat, zmat3));

}
