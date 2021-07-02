#pragma once

#include <algorithm>

#include <lila/vector.h>
#include <lila/matrix.h>
#include <lila/views/vector_view.h>
#include <lila/views/matrix_view.h>

namespace lila {

template <class coeff_t> class Vector;
template <class coeff_t> class VectorView;
template <class coeff_t> class Matrix;
template <class coeff_t> class MatrixView;

template <class coeff_t, class function_t>
inline void Map(Vector<coeff_t> &v, function_t func) {
  std::for_each(v.begin(), v.end(), func);
}

template <class coeff_t, class function_t>
inline void Map(VectorView<coeff_t> v, function_t func) {
  coeff_t *data = v.data();
  size_type inc = v.inc();
  size_type size = v.n() * inc;

  for (size_type i = 0; i < size; i += inc)
    func(data[i]);
}

template <class coeff_t, class function_t>
inline void Map(Matrix<coeff_t> &A, function_t func) {
  std::for_each(A.begin(), A.end(), func);
}

template <class coeff_t, class function_t>
inline void Map(MatrixView<coeff_t> A, function_t func) {
  size_type ld = A.ld();
  size_type incm = A.incm();
  size_type incn = A.incn();
  size_type sizem = A.m() * incm;
  size_type sizen = A.n() * incn;

  for (size_type idxn = 0; idxn < sizen; idxn += incn) {
    coeff_t *data = A.data() + idxn * ld;
    for (size_type idxm = 0; idxm < sizem; idxm += incm)
      func(data[idxm]);
  }
}

} // namespace lila
