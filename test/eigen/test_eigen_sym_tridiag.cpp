#include "../catch.hpp"
#include <lila/all.h>

template <class coeff_t> void test_eigen_sym_tridiag() {
  using namespace lila;

  int n = 10;

  // Test Eigen
  for (int seed : lila::range<int>(10)) {

    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, seed);

    Vector<coeff_t> diag(n);
    Vector<coeff_t> offdiag(n - 1);
    Random(diag, fgen);
    Random(offdiag, fgen);

    auto fullmat = Matrix<coeff_t>(n, n);
    for (int i = 0; i < n; ++i)
      fullmat(i, i) = diag(i);
    for (int i = 0; i < n - 1; ++i) {
      fullmat(i, i + 1) = offdiag(i);
      fullmat(i + 1, i) = offdiag(i);
    }

    auto [eigs, evecs] = EigenSymTridiag(diag, offdiag);
    auto eigs2 = EigenvaluesSymTridiag(diag, offdiag);
    REQUIRE(close(eigs, eigs2));

    // Check if indeed eigenvalues
    for (int i=0; i<evecs.ncols(); ++i) {
      Vector<coeff_t> evec = evecs(ALL, i);
      coeff_t test = Dot(evec, Mult(fullmat, evec));
      REQUIRE(close(eigs(i), test));
    }
  }
}

TEST_CASE("eigen_sym_tridiag", "[decomp]") {
  lila::Log("Test eigen_sym_tridiag");

  test_eigen_sym_tridiag<float>();
  test_eigen_sym_tridiag<double>();
}
