#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t>
lila::Matrix<coeff_t> simple_multiply(const lila::Matrix<coeff_t> &A,
                                      const lila::Matrix<coeff_t> &B) {
  lila::Matrix<coeff_t> C(A.nrows(), B.ncols());
  lila::Zeros(C);
  for (int i = 0; i < C.nrows(); ++i)
    for (int j = 0; j < C.ncols(); ++j)
      for (int k = 0; k < A.ncols(); ++k)
        C(i, j) += A(i, k) * B(k, j);
  return C;
}

template <class coeff_t>
lila::Vector<coeff_t> simple_multiply(const lila::Matrix<coeff_t> &A,
                                      const lila::Vector<coeff_t> &X) {
  lila::Vector<coeff_t> Y(A.nrows());
  lila::Zeros(Y);
  for (int i = 0; i < A.nrows(); ++i)
    for (int j = 0; j < A.ncols(); ++j)
      Y(i) += A(i, j) * X(j);
  return Y;
}

template <class coeff_t> void test_mult() {
  int m = 7;
  int k = 9;
  int n = 11;

  int seed = 42;

  lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
  lila::uniform_gen_t<coeff_t> fgen(fdist, seed);

  lila::Matrix<coeff_t> A(m, k);
  lila::Matrix<coeff_t> B(k, n);
  lila::Matrix<coeff_t> C(m, n);
  lila::Vector<coeff_t> X(k);
  lila::Vector<coeff_t> Y(m);
  Random(A, fgen);
  Random(B, fgen);
  Random(C, fgen);
  Random(X, fgen);
  Mult(A, B, C);

  lila::Matrix<coeff_t> C2 = simple_multiply(A, B);
  REQUIRE(lila::close<coeff_t>(C, C2));
  Mult(A, X, Y);
  lila::Vector<coeff_t> Y2 = simple_multiply(A, X);

  REQUIRE(lila::close<coeff_t>(Y, Y2));

  Mult(B, A, C, coeff_t(.5), coeff_t(.6), 'T', 'T');

  A.resize(m, k);
  B.resize(n, k);
  Mult(A, B, C, coeff_t(.5), coeff_t(.6), 'N', 'T');

  A.resize(k, m);
  B.resize(k, n);
  Mult(A, B, C, coeff_t(.5), coeff_t(.6), 'T', 'N');

  A.resize(2, 3);
  Zeros(A);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(1, 0) = 3;
  A(1, 1) = 4;
  B.resize(2, 2);
  B(0, 0) = 0;
  B(0, 1) = 5;
  B(1, 0) = 6;
  B(1, 1) = 7;
  Kron(A, B, C);
  // Print(C);
}

TEST_CASE("mult", "[algebra]") {
  lila::Log("Test mult");

  test_mult<float>();
  test_mult<double>();
  test_mult<std::complex<float>>();
  test_mult<std::complex<double>>();
}
