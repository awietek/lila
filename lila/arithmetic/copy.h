#pragma once

#include <cassert>

#include <lila/blaslapack/blaslapack.h>
#include <lila/matrix.h>
#include <lila/vector.h>
#include <lila/views/matrix_view.h>
#include <lila/views/vector_view.h>

namespace lila {

template <class coeff_t> class Vector;
template <class coeff_t> class Matrix;
template <class coeff_t> class VectorView;
template <class coeff_t> class MatrixView;

template <class coeff_t>
inline void Copy(VectorView<coeff_t> &&X, VectorView<coeff_t> &&Y) {
  using size_type = blaslapack::blas_size_t;

  assert(X.N() == Y.N()); // Check if valid dimensions
  size_type N = X.N();
  size_type incx = X.inc();
  size_type incy = Y.inc();
  blaslapack::copy(&N, LILA_BLAS_CONST_CAST(coeff_t, X.data()), &incx,
                   LILA_BLAS_CAST(coeff_t, Y.data()), &incy);
}

template <class coeff_t>
inline void Copy(Matrix<coeff_t> const &X, Matrix<coeff_t> &Y) {
  using size_type = blaslapack::blas_size_t;
  size_type dx = X.size();
  size_type dy = Y.size();
  assert(dx == dy); // Check if valid dimensions
  size_type inc = 1;
  blaslapack::copy(&dx, LILA_BLAS_CONST_CAST(coeff_t, X.data()), &inc,
                   LILA_BLAS_CAST(coeff_t, Y.data()), &inc);
}

template <class coeff_t>
inline void Copy(MatrixView<coeff_t> && X, MatrixView<coeff_t> && Y) {
  using size_type = blaslapack::blas_size_t;
  std::cout << X.M() << " " << Y.M() << "\n";
  assert(X.M() == Y.M());
  assert(X.N() == Y.N());
  size_type M = X.M();
  size_type N = X.N();
  size_type ldX = X.ld();
  size_type ldY = Y.ld();
  size_type inc = 1;

  // Perform a column-wise copy
  for (int col = 0; col < N; ++col) {
    blaslapack::copy(&M, LILA_BLAS_CONST_CAST(coeff_t, X.data() + col * ldX),
                     &inc, LILA_BLAS_CAST(coeff_t, Y.data() + col * ldY), &inc);
  }
}

} // namespace lila
