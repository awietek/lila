#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t> void test_add() {
  using namespace lila;
  int m = 7;
  int n = 11;
  int step = 2;

  // Vector
  auto v1 = Random<coeff_t>(n);
  auto v2 = ZerosLike(v1);
  auto v3 = ZerosLike(v1);
  v2 += v1;
  v3 = v2 + v1;
  for (int i = 0; i < n; ++i) {
    REQUIRE(close(v1(i), v2(i)));
    REQUIRE(close(v1(i), v3(i) / (coeff_t)2.));
  }

  // VectorView (contiguous)
  v1 = Random<coeff_t>(n);
  v2 = ZerosLike(v1);
  v3 = ZerosLike(v1);
  v2({1, n - 1}) += v1({1, n - 1});
  v3({1, n - 1}) = v2({1, n - 1}) + v1({1, n - 1});
  for (int i = 1; i < n - 1; ++i) {
    REQUIRE(close(v1(i), v2(i)));
    REQUIRE(close(v1(i), v3(i) / (coeff_t)2.));
  }

  // VectorView (steps)
  v1 = Random<coeff_t>(n);
  v2 = ZerosLike(v1);
  v3 = ZerosLike(v1);
  v2({1, n - 1, step}) += v1({1, n - 1, step});
  v3({1, n - 1, step}) = v2({1, n - 1, step}) + v1({1, n - 1, step});
  for (int i = 1; i < n - 1; i += step) {
    REQUIRE(close(v1(i), v2(i)));
    REQUIRE(close(v1(i), v3(i) / (coeff_t)2.));
  }

  // Matrix
  auto A1 = Random<coeff_t>(m, n);
  auto A2 = ZerosLike(A1);
  auto A3 = ZerosLike(A1);
  A2 += A1;
  A3 = A2 + A1;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j) {
      REQUIRE(close(A1(i, j), A2(i, j)));
      REQUIRE(close(A1(i, j), A3(i, j) / (coeff_t)2.));
    }

  // MatrixView (contiguous)
  A1 = Random<coeff_t>(m, n);
  A2 = ZerosLike(A1);
  A2({1, n - 1}, {1, n - 1}) += A1({1, n - 1}, {1, n - 1});
  A3({1, n - 1}, {1, n - 1}) =
      A2({1, n - 1}, {1, n - 1}) + A1({1, n - 1}, {1, n - 1});
  for (int i = 1; i < m - 1; ++i)
    for (int j = 1; j < n - 1; ++j) {
      REQUIRE(close(A1(i, j), A2(i, j)));
      REQUIRE(close(A1(i, j), A3(i, j) / (coeff_t)2.));
    }
  // MatrixView (steps)
  A1 = Random<coeff_t>(m, n);
  A2 = ZerosLike(A1);
  A2({1, n - 1, step}, {1, n - 1, step}) +=
      A1({1, n - 1, step}, {1, n - 1, step});
  A3({1, n - 1, step}, {1, n - 1, step}) =
      A2({1, n - 1, step}, {1, n - 1, step}) +
      A1({1, n - 1, step}, {1, n - 1, step});
  for (int i = 1; i < m - 1; i += step)
    for (int j = 1; j < n - 1; j += step) {
      REQUIRE(close(A1(i, j), A2(i, j)));
      REQUIRE(close(A1(i, j), A3(i, j) / (coeff_t)2.));
    }
}

TEST_CASE("add", "[arithmetic]") {
  test_add<float>();
  test_add<double>();
  test_add<std::complex<float>>();
  test_add<std::complex<double>>(); 
}
