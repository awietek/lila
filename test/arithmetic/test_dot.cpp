#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t>
void test_dot()
{
  using namespace lila;
  int n=20;
  
  // Test Dot
  for (int seed : lila::range<int>(10)) {
  
    lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
    lila::uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    lila::Vector<coeff_t> v1(n);
    Random(v1, fgen);
   
    lila::Vector<coeff_t> v2(n);
    Random(v2, fgen);

    coeff_t dot_lila = Dot(v1, v2);
    coeff_t dot_test = 0;
    for (int i = 0; i < n; ++i)
      dot_test += lila::conj(v1(i)) * v2(i);

    REQUIRE(close(dot_lila, dot_test));
  }
}


TEST_CASE( "dot", "[arithmetic]" ) {
  test_dot<float>();
  test_dot<double>();
  test_dot<std::complex<float>>();
  test_dot<std::complex<double>>();
}
