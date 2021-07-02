#pragma once

namespace lila {
template <class coeff_t> coeff_t Trace(Matrix<coeff_t> const &mat) {
  // TODO: optimize ??
  coeff_t res = 0;
  for (size_type i = 0; i < std::min(mat.nrows(), mat.ncols()); ++i)
    res += mat(i, i);
  return res;
}

template <class coeff_t> Matrix<coeff_t> Transpose(Matrix<coeff_t> const &mat) {
  Matrix<coeff_t> mat_t(mat.ncols(), mat.nrows());
  for (int i = 0; i < mat.nrows(); ++i)
    for (int j = 0; j < mat.ncols(); ++j)
      mat_t(j, i) = mat(i, j);
  return mat_t;
}

template <class coeff_t> Vector<coeff_t> Diag(Matrix<coeff_t> const &mat) {
  long size = std::min(mat.nrows(), mat.ncols());
  Vector<coeff_t> vec_t(size);
  for (long i = 0; i < size; ++i)
    vec_t(i) = mat(i, i);
  return vec_t;
}

} // namespace lila
