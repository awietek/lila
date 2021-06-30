#pragma once

#include <cassert>

#include <lila/arithmetic/copy.h>
#include <lila/views/matrix_view.h>
#include <lila/views/vector_view.h>

namespace lila {

template <class coeff_t, class coeff_t2>
inline void Assign(VectorView<coeff_t> const &X, coeff_t2 val2) {
  coeff_t val = (coeff_t)val2;
  coeff_t *ptr = X.data();
  for (size_type i = 0; i < X.N(); ++i) {
    *ptr = val;
    ptr += X.inc();
  }
}

template <class coeff_t>
inline void Assign(VectorView<coeff_t> const &X, VectorView<coeff_t> const &Y) {
  Copy(Y, X);
}

template <class coeff_t>
inline void Assign(VectorView<coeff_t> const &X, Vector<coeff_t> const &Y) {
  Assign(X, VectorView<coeff_t>(const_cast<Vector<coeff_t>&>(Y)));
}


template <class coeff_t, class coeff_t2>
inline void Assign(MatrixView<coeff_t> const &X, coeff_t2 val2) {
  coeff_t val = (coeff_t)val2;
  coeff_t *ptr = X.data();
  for (size_type j = 0; j < X.N(); ++j) {
    for (size_type i = 0; i < X.M(); ++i) {
      *ptr = val;
      ++ptr;
    }
    ptr += X.ld() - X.M();
  }
}

template <class coeff_t>
inline void Assign(MatrixView<coeff_t> const &X, MatrixView<coeff_t> const &Y) {
  Copy(Y, X);
}

template <class coeff_t>
inline void Assign(MatrixView<coeff_t> const &X, Matrix<coeff_t> const &Y) {
  Assign(X, MatrixView<coeff_t>(const_cast<Matrix<coeff_t>&>(Y)));
}


template <class coeff_t, class coeff_t2>
VectorView<coeff_t> &operator^=(VectorView<coeff_t> &&X, coeff_t2 val) {
  Assign(X, val);
  return X;
}

template <class coeff_t>
VectorView<coeff_t> &operator^=(VectorView<coeff_t> &&X,
                                VectorView<coeff_t> const &Y) {
  Assign(X, Y);
  return X;
}

template <class coeff_t, class coeff_t2>
MatrixView<coeff_t> &operator^=(MatrixView<coeff_t> &&X, coeff_t2 val) {
  Assign(X, val);
  return X;
}

template <class coeff_t>
MatrixView<coeff_t> &operator^=(MatrixView<coeff_t> &&X,
                                MatrixView<coeff_t> const &Y) {
  Assign(X, Y);
  return X;
}

} // namespace lila
