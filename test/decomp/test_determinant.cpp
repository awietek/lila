#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t> void test_determinant(int m) {

  // Test direct solve
  for (int seed : lila::range<int>(10)) {

    lila::uniform_dist_t<coeff_t> dist(-1., 1.);
    lila::uniform_gen_t<coeff_t> gen(dist, seed);

    auto A = Random(m, m, gen);
    lila::complex_t<coeff_t> determinant = Determinant(A);
    auto eigenvalues = Eigenvalues(A);
    lila::complex_t<coeff_t> prod_eigs = 1.0;
    for (auto eig : eigenvalues)
      prod_eigs *= eig;
    REQUIRE(lila::close(determinant, prod_eigs));
  }
}

TEST_CASE("determinant", "[decomp]") {
  lila::Log("Test determinant");

  for (int m = 1; m < 8; ++m) {
    test_determinant<float>(m);
    test_determinant<double>(m);
    test_determinant<std::complex<float>>(m);
    test_determinant<std::complex<double>>(m);
  }
}
