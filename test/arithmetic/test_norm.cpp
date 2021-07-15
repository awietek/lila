#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t>
void test_norm()
{
  using namespace lila;
  int n=20;
  
  // Test Norm
  for (int seed : lila::range<int>(10)) {
  
    lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
    lila::uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    lila::Vector<coeff_t> v1(n);
    Random(v1, fgen);
   
    coeff_t norm_lila = Norm(v1);
    coeff_t norm_test = 0;
    for (int i = 0; i < n; ++i)
      norm_test += lila::conj(v1(i)) * v1(i);

    LilaPrint(norm_lila);
    LilaPrint(std::sqrt(norm_test));
    REQUIRE(close(norm_lila, std::sqrt(norm_test)));
  }
}


TEST_CASE( "norm", "[arithmetic]" ) {
  test_norm<float>();
  test_norm<double>();
  test_norm<std::complex<float>>();
  test_norm<std::complex<double>>();
}
