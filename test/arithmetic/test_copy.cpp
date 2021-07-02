#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t> void test_copy() {
  using namespace lila;
  int m = 7;
  int n = 11;
  int step = 2;

  // Vector
  auto v1 = Random<coeff_t>(n);
  auto v2 = ZerosLike(v1);
  Copy(v1, v2);
  for (int i = 0; i < n; ++i)
    REQUIRE(v1(i) == v2(i));

  // VectorView (contiguous)
  v1 = Random<coeff_t>(n);
  v2 = ZerosLike(v1);
  Copy(v1({1, n - 1}), v2({1, n - 1}));
  LilaPrint(v1);
  LilaPrint(v2);
  for (int i = 1; i < n - 1; ++i)
    REQUIRE(v1(i) == v2(i));

  // VectorView (steps)
  v1 = Random<coeff_t>(n);
  v2 = ZerosLike(v1);
  Copy(v1({1, n - 1, step}), v2({1, n - 1, step}));

  LilaPrint(v1);
  LilaPrint(v2);
  for (int i = 1; i < n - 1; i += step)
    REQUIRE(v1(i) == v2(i));

  // Matrix
  auto A1 = Random<coeff_t>(m, n);
  auto A2 = ZerosLike(A1);
  Copy(A1, A2);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      REQUIRE(A1(i, j) == A2(i, j));

  // MatrixView (contiguous)
  A1 = Random<coeff_t>(m, n);
  A2 = ZerosLike(A1);
  Copy(A1({1, m - 1}, {1, n - 1}), A2({1, m - 1}, {1, n - 1}));

  LilaPrint(A1);
  LilaPrint(A2);
  for (int i = 1; i < m - 1; ++i)
    for (int j = 1; j < n - 1; ++j)
      REQUIRE(A1(i, j) == A2(i, j));

  // MatrixView (steps)
  A1 = Random<coeff_t>(m, n);
  A2 = ZerosLike(A1);
  Copy(A1({1, m - 1, step}, {1, n - 1, step}),
       A2({1, m - 1, step}, {1, n - 1, step}));
  for (int i = 1; i < m - 1; i += step)
    for (int j = 1; j < n - 1; j += step)
      REQUIRE(A1(i, j) == A2(i, j));
}

TEST_CASE("copy", "[arithmetic]") {
  test_copy<float>();
  test_copy<double>();
  test_copy<std::complex<float>>();
  test_copy<std::complex<double>>();
}
