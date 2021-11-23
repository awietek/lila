#pragma once

#include <algorithm>

#include <lila/matrix.h>
#include <lila/vector.h>
#include <lila/views/matrix_view.h>
#include <lila/views/vector_view.h>

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
  lila_size_t inc = v.inc();
  lila_size_t size = v.n() * inc;

  for (lila_size_t i = 0; i < size; i += inc)
    func(data[i]);
}

template <class coeff_t, class function_t>
inline void Map(Matrix<coeff_t> &A, function_t func) {
  std::for_each(A.begin(), A.end(), func);
}

template <class coeff_t, class function_t>
inline void Map(MatrixView<coeff_t> A, function_t func) {
  lila_size_t ld = A.ld();
  lila_size_t incm = A.incm();
  lila_size_t incn = A.incn();
  lila_size_t sizem = A.m() * incm;
  lila_size_t sizen = A.n() * incn;

  for (lila_size_t idxn = 0; idxn < sizen; idxn += incn) {
    coeff_t *data = A.data() + idxn * ld;
    for (lila_size_t idxm = 0; idxm < sizem; idxm += incm)
      func(data[idxm]);
  }
}

} // namespace lila
