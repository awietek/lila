#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t>
void test_solve()
{
  int n=5;
  int k=3;
  
  // Test direct solve
  for (int seed : lila::range<int>(10)) {
  
    lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
    lila::uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    lila::Matrix<coeff_t> A(n, n);
    lila::Matrix<coeff_t> B(n, k);
    lila::Vector<coeff_t> b(n);

    Random(A, fgen);
    Random(B, fgen);
    Random(b, fgen);

    // Solve with multiple right hand sides
    lila::Matrix<coeff_t> A_save(A);
    lila::Matrix<coeff_t> X(B);
    Solve(A, X);
    lila::Matrix<coeff_t> Prod(n, k);
    Mult(A_save, X, Prod);
    REQUIRE(lila::close<coeff_t>(B, Prod));
    // Print(B);
    // Print(Prod);
 
    // Solve using explicit LU Decomposition
    A = A_save;
    auto ipiv = LUDecompose(A);
    X = B;
    LUSolve(A, ipiv, X);
    Zeros(Prod);
    Mult(A_save, X, Prod);
    REQUIRE(lila::close<coeff_t>(B, Prod));

    // Solve for one right hand side
    A = A_save;
    lila::Vector<coeff_t> x(b);
    Solve(A, x);
    lila::Vector<coeff_t> prod(n);
    Mult(A_save, x, prod);
    REQUIRE(lila::close<coeff_t>(b, prod));

    // Solve one rhs explicit LU Decomposition
    A = A_save;
    ipiv = LUDecompose(A);
    x = b;
    LUSolve(A, ipiv, x);
    Zeros(prod);
    Mult(A_save, x, prod);
    REQUIRE(lila::close<coeff_t>(B, Prod));

    // Test inversion
    A = A_save;
    lila::Matrix<coeff_t> A_inv = A;
    Invert(A_inv);
    lila::Matrix<coeff_t> Id(A);
    Identity(Id);

    Zeros(Prod);
    Mult(A_inv, A, Prod);
    REQUIRE(lila::close<coeff_t>(Id, Prod));

    Zeros(Prod);
    Mult(A, A_inv, Prod);
    REQUIRE(lila::close<coeff_t>(Id, Prod));
  }
  
}


TEST_CASE( "solve", "[decomp]" ) {
  test_solve<float>();
  test_solve<double>();
  test_solve<std::complex<float>>();
  test_solve<std::complex<double>>();
}
