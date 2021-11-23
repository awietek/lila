#pragma once

#include <math.h>

#include <lila/matrix.h>
#include <lila/arithmetic/add.h>
#include <lila/algebra/mult.h>
#include <lila/blaslapack/blaslapack.h>
#include <lila/eigen/eigen_sym.h>
#include <lila/numeric/complex.h>
#include <lila/special/special.h>

namespace lila {

template <class coeff_t, class function_t>
inline void FunctionSym(Matrix<coeff_t> &matrix, function_t fun,
                        char uplo = 'U') {
  assert(matrix.m() == matrix.n());

  auto Qmat = matrix;
  auto eigs = EigenSymInplace(Qmat, true, uplo);
  auto Bmat = Qmat;

  // Multiply with diagonal matrix, where function has been applied
  for (lila_size_t j = 0; j < eigs.n(); ++j) {
    coeff_t fun_of_eig = static_cast<coeff_t>(eigs(j));
    fun(fun_of_eig);

    // Scale j-th column
    // blaslapack::scal(&n, LILA_BLAS_CAST(coeff_t, &fun_of_eig),
    //                  LILA_BLAS_CAST(coeff_t, Bmat.data()) + j * n, &incx);
    Bmat(ALL, j) *= fun_of_eig;
  }
  Zeros(matrix);
  Mult(Bmat, Qmat, matrix, (coeff_t)1., (coeff_t)0., 'N', 'C');
}

template <class coeff_t>
inline void ExpSym(Matrix<coeff_t> &matrix, coeff_t alpha = 1.,
                   char uplo = 'U') {
  FunctionSym(
      matrix, [alpha](coeff_t &x) { x = exp(alpha * x); }, uplo);
}

template <class coeff_t>
inline void LogSym(Matrix<coeff_t> &matrix, coeff_t alpha = 1.,
                   char uplo = 'U') {
  FunctionSym(
      matrix, [alpha](coeff_t &x) { x = log(x) / alpha; }, uplo);
}

} // namespace lila
