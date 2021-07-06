#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t>
void test_expm()
{
  using namespace lila;
  int n=10;
  coeff_t alpha=1.23;
  for (int seed : range<int>(10)) {
  
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    Matrix<coeff_t> A(n, n);
    Random(A, fgen);
    // LilaPrint(A);

    auto Aexp = ExpM(A, alpha);
    // LilaPrint(Aexp);


    // Check whether hermitian gives same result
    auto A1 = A;
    Add(Herm(A1), A1);
    // LilaPrint(A1);

    // auto LILA_CLK(t1);
    Aexp = ExpM(A1, alpha);
    // auto LILA_CLK(t2);
    // LILA_TIME_MICRO("pade", t1, t2);
    // LilaPrint(Aexp);


    auto Aexp2 = A1;
    // LILA_CLK(t1);
    ExpSym(Aexp2, alpha);
    // LILA_CLK(t2);
    // LILA_TIME_MICRO("herm", t1, t2);
    // LilaPrint(Aexp);
    // LilaPrint(Aexp2);
    REQUIRE(close(Aexp2, Aexp));

    // Pade from EXPOKIT is actually really good in timing!!
    
  }
}


TEST_CASE( "expm", "[functions]" ) {
  test_expm<float>();
  test_expm<double>();
  test_expm<std::complex<float>>();
  test_expm<std::complex<double>>();
}
