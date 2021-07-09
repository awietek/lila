#pragma once

#include <lila/decomp/solve.h>
#include <lila/matrix.h>

#include <lila/utils/print.h>

namespace lila {
template <class coeff_t> coeff_t DeterminantInplace(Matrix<coeff_t> &A) {
  assert(A.m() == A.n());
  auto ipiv = LUDecompose(A);

  coeff_t det = 1.0;
  int sign = 0;
  for (int i = 0; i < A.m(); ++i) {

    if (ipiv[i] != i + 1) {
      sign ^= 1;
    }

    det *= A(i, i);
  }
  if (sign & 1)
    return -det;
  else
    return det;
}

template <class coeff_t> coeff_t Determinant(Matrix<coeff_t> A) {
  return DeterminantInplace(A);
}

} // namespace lila
