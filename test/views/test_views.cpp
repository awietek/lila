#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t>
void test_views()
{
  using namespace lila;

  auto vec = Random<coeff_t>(8);
  auto subvec = vec({2,8,2});
  LilaPrint(vec);
  LilaPrint(subvec);

  auto mat = Random<coeff_t>(5,5);
  auto submat = mat({2,4}, {2,4});
  LilaPrint(mat);
  LilaPrint(submat);
}


TEST_CASE( "views", "[core]" ) {
  test_views<float>();
  test_views<double>();
  test_views<std::complex<float>>();
  test_views<std::complex<double>>();
}
