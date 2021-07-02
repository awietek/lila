#include "catch.hpp"

#include <lila/all.h>

template <class coeff_t>
void test_resize()
{
  using namespace lila;

  int m=6;
  int n=5;

  // Test Resize
  for (int m_new=1; m_new < 10; ++m_new)
    for (int n_new=1; n_new < 10; ++n_new)
      {
	for (int seed : lila::range<int>(3)) 
	  {
      
	    uniform_dist_t<coeff_t> dist(-1., 1.);
	    uniform_gen_t<coeff_t> gen(dist, seed);
      
	    Matrix<coeff_t> A(m, n);
	    Random(A, gen);
	    auto A_res = A;
	    A_res.resize(m_new, n_new);
	    for (int i=0; i < m_new; ++i)
	      for (int j=0; j < n_new; ++j)
		{
		  if ((i < std::min(m_new, m)) && (j < std::min(n_new, n))) 
		    REQUIRE(A_res(i, j) == A(i, j));
		  else REQUIRE(close<coeff_t>(A_res(i, j), 0.));
		}
      	  }
      }
}


TEST_CASE( "resize", "[core]" ) {
  test_resize<float>();
  test_resize<double>();
  test_resize<std::complex<float>>();
  test_resize<std::complex<double>>();
}
