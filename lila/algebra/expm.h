#pragma once

#include <lila/matrix.h>

namespace lila {

template <class coeff_t>
inline real_t<coeff_t> NormLi(Matrix<coeff_t> const &X) {
  real_t<coeff_t> value = 0.0;
  for (int i = 0; i < X.nrows(); i++) {
    real_t<coeff_t> row_sum = 0.0;
    for (int j = 0; j < X.ncols(); j++) {
      row_sum += std::abs(X(i, j));
    }
    value = std::max(value, row_sum);
  }
  return value;
}

template <class coeff_t> inline real_t<coeff_t> log2abs(coeff_t x) {
  real_t<coeff_t> value;

  if (x == 0.0) {
    value = -1e30;
  } else {
    value = log(std::abs(x)) / log(2.0);
  }

  return value;
}

template <class coeff_t>
inline Matrix<coeff_t> ExpM(Matrix<coeff_t> const &A, coeff_t alpha = 1.) {
  int n = A.nrows();
  assert(n == A.ncols());
  const int q = 6;
  Matrix<coeff_t> a2 = alpha * A;
  real_t<coeff_t> a_norm = NormLi(a2);
  int ee = (int)log2abs(a_norm) + 1;
  int s = std::max(0, ee + 1);
  coeff_t t = 1.0 / pow(2.0, s);
  Scale(t, a2);
  Matrix<coeff_t> x = a2;
  coeff_t c = 0.5;
  Matrix<coeff_t> e = Identity<coeff_t>(n);
  Add(a2, e, c);
  Matrix<coeff_t> d = Identity<coeff_t>(n);
  Add(a2, d, -c);
  int p = 1;
  for (int k = 2; k <= q; k++) {
    c = c * (coeff_t)(q - k + 1) / (coeff_t)(k * (2 * q - k + 1));
    x = Mult(a2, x);  //SEGFAULT WHY?
    Add(x, e, c);
    if (p) {
      Add(x, d, c);
    } else {
      Add(x, d, -c);
    }
    p = !p;
  }

  Solve(d, e);
  for (int k = 1; k <= s; k++) {
    e = Mult(e, e);
  }
  return e;
}

} // namespace lila
