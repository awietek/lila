#pragma once

#include <cassert>
#include <lila/blaslapack/blaslapack.h>

namespace lila {

template <class coeff_t>
inline coeff_t Dot(Vector<coeff_t> const &v, Vector<coeff_t> const &w) {
  assert(v.n() == w.n());
  blas_size_t n = v.n();
  blas_size_t inc = 1;
  blaslapack::blas_t<coeff_t> result =
      blaslapack::dot(&n, LILA_BLAS_CONST_CAST(coeff_t, v.data()), &inc,
                      LILA_BLAS_CONST_CAST(coeff_t, w.data()), &inc);
  return blaslapack::blas_to_lila(result);
}

template <class coeff_t>
inline coeff_t Dot(VectorView<coeff_t> const &v, VectorView<coeff_t> const &w) {
  assert(v.n() == w.n());
  blas_size_t n = v.n();
  blas_size_t v_inc = v.inc();
  blas_size_t w_inc = w.inc();
  blaslapack::blas_t<coeff_t> result =
      blaslapack::dot(&n, LILA_BLAS_CONST_CAST(coeff_t, v.data()), &v_inc,
                      LILA_BLAS_CONST_CAST(coeff_t, w.data()), &w_inc);
  return blaslapack::blas_to_lila(result);
}

template <class coeff_t>
inline coeff_t Dot(Vector<coeff_t> const &v, VectorView<coeff_t> const &w) {
  return Dot(VectorView<coeff_t>(v), w);
}

template <class coeff_t>
inline coeff_t Dot(VectorView<coeff_t> const &v, Vector<coeff_t> const &w) {
  return Dot(v, VectorView<coeff_t>(w));
}

} // namespace lila
