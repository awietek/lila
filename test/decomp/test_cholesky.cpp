#include "../catch.hpp"

#include <lila/all.h>

#include <iostream>
template <class coeff_t>
void test_cholesky()
{
  using namespace lila;
  int n=10;
  
  // Test Cholesky
  for (int seed=0; seed<10; ++seed) 
    {
      Matrix<coeff_t> A(n, n);
    
      // Create positive definite Matrix
      std::vector<lila::Vector<coeff_t>> vecs; 
      for (int k=0; k<n; ++k)
	vecs.push_back(Random<coeff_t>(n));
      for (int i=0; i<n; ++i)
	for (int j=0; j<n; ++j)
	  A(i,j) = Dot(vecs[i], vecs[j]);
    
      // LilaPrint(A);
      auto U = Cholesky(A, 'U');
      // LilaPrint(U);
      auto produ = Mult(Herm(U), U);
      // LilaPrint(produ);
      REQUIRE(close(produ, A));      

      auto L = Cholesky(A, 'L');
      // LilaPrint(L);
      auto prodl = Mult(L, Herm(L));
      // LilaPrint(prodl);
      REQUIRE(close(prodl, A)); 

    }
}


TEST_CASE( "cholesky", "[decomp]" ) {
  test_cholesky<float>();
  test_cholesky<double>();
  test_cholesky<std::complex<float>>();
  test_cholesky<std::complex<double>>();
}
