#pragma once

#include <utility>

#include <lila/blaslapack/blaslapack.h>
#include <lila/vector.h>

namespace lila {

template <class coeff_t>
inline void EigenvaluesSymTridiagInplace(Vector<coeff_t> &diag,
                                         Vector<coeff_t> &offdiag) {
  assert(diag.size() <= offdiag.size() + 1);
  char job = 'N';
  blas_size_t N = diag.size();
  blas_size_t ldz = 1;
  blas_size_t info = 0;
  blaslapack::stev(&job, &N, diag.data(), offdiag.data(), NULL, &ldz, NULL,
                   &info);
}

template <class coeff_t>
inline Vector<coeff_t> EigenvaluesSymTridiag(Vector<coeff_t> const &diag,
                                             Vector<coeff_t> const &offdiag) {
  auto diag_copy = diag;
  auto offdiag_copy = offdiag;
  EigenvaluesSymTridiagInplace(diag_copy, offdiag_copy);
  return diag_copy;
}

template <class coeff_t>
inline Matrix<coeff_t> EigenSymTridiagInplace(Vector<coeff_t> &diag,
                                              Vector<coeff_t> &offdiag) {
  assert(diag.size() <= offdiag.size() + 1);

  char job = 'V';
  blas_size_t N = diag.size();
  blas_size_t ldz = N;
  blas_size_t info = 0;

  auto eigenvectors = Matrix<coeff_t>(N, N);

  Vector<coeff_t> work(2 * N - 2);
  blaslapack::stev(&job, &N, diag.data(), offdiag.data(), eigenvectors.data(),
                   &ldz, work.data(), &info);
  return eigenvectors;
}

template <class coeff_t>
inline std::pair<Vector<coeff_t>, Matrix<coeff_t>>
EigenSymTridiag(Vector<coeff_t> const &diag, Vector<coeff_t> const &offdiag) {
  auto eigenvalues = diag;
  auto offdiag_copy = offdiag;
  auto eigenvectors = EigenSymTridiagInplace(eigenvalues, offdiag_copy);
  return {eigenvalues, eigenvectors};
}

} // namespace lila
