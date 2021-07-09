#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t>
void test_matrixfunction()
{
  using namespace lila;
  int n=5;
  coeff_t alpha=1.23;
  for (int seed : range<int>(10)) {
  
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    Matrix<coeff_t> A(n, n);
    Random(A, fgen);
    auto A1 = A;
    Add(Herm(A1), A1);
    auto Aexp = A1;
    ExpSym(Aexp, alpha);
    auto Aret = Aexp;
    LogSym(Aret, alpha);

    // LilaPrint(A1);
    // LilaPrint(Aret);
    REQUIRE(close(A1, Aret));

    auto v2 = EigenvaluesSym(Aexp);
    for (auto e : v2)
      REQUIRE(e >= 0);
    
  }
}


TEST_CASE( "matrixfunction", "[function]" ) {
  // test_matrixfunction<float>();
  test_matrixfunction<double>();
  // test_matrixfunction<std::complex<float>>();
  test_matrixfunction<std::complex<double>>();
}
