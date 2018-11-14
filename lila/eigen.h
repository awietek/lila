// Copyright 2018 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LILA_EIGEN_H_
#define LILA_EIGEN_H_

#include "complex.h"

namespace lila {

  template <class coeff_t>
  struct EigenResults {
    Vector<complex_t<coeff_t>> eigenvalues;
    Matrix<coeff_t> left_eigenvectors;
    Matrix<coeff_t> right_eigenvectors;    
  };

  template <class coeff_t>
  inline EigenResults<coeff_t> 
  Eigen(Matrix<coeff_t>& A, bool do_right_eigenvectors=false, 
	bool do_left_eigenvectors=false)
  {
    using size_type = blaslapack::blas_size_t;
    assert(A.nrows() == A.ncols());

    const char jobvl = do_left_eigenvectors ? 'V' : 'N';
    const char jobvr = do_right_eigenvectors ? 'V' : 'N';
    const size_type n = A.nrows();
    const size_type lda = n;
 
    const size_type ldvl = do_left_eigenvectors ? n : 1;
    const size_type ldvr = do_right_eigenvectors ? n : 1;
    
    EigenResults<coeff_t> result;
    result.eigenvalues.resize(n);
    if (do_left_eigenvectors) result.left_eigenvectors.resize(ldvl, n);
    if (do_right_eigenvectors) result.right_eigenvectors.resize(ldvr, n);

    // get optimal work size
    Vector<complex_t<coeff_t>>& w = result.eigenvalues;
    size_type lwork = -1;
    std::vector<coeff_t> work(1);
    int info = 0;
    blaslapack::geev(&jobvl, &jobvr, &n, A.data(), &lda, w.data(), 
		     result.left_eigenvectors.data(), &ldvl,
		     result.right_eigenvectors.data(), &ldvr, 
		     work.data(), &lwork, &info);
    assert(info == 0);
    lwork = static_cast<size_type>(real(work[0]));
    work.resize(lwork);

    // Run eigenvalue computation
    blaslapack::geev(&jobvl, &jobvr, &n, A.data(), &lda, w.data(), 
		     result.left_eigenvectors.data(), &ldvl,
		     result.right_eigenvectors.data(), &ldvr, 
		     work.data(), &lwork, &info);
    assert(info == 0);
    return result;
  }

  template <class coeff_t>
  inline Vector<complex_t<coeff_t>> 
  Eigenvalues(Matrix<coeff_t>& A) {
    EigenResults<coeff_t> res = Eigen(A);
    return res.eigenvalues;
  }

  template <class coeff_t>
  inline Vector<real_t<coeff_t>> 
  EigenH(Matrix<coeff_t>& A, bool do_eigenvectors=true, char uplo='U')
  {
    using size_type = blaslapack::blas_size_t;
    assert(A.nrows() == A.ncols());

    const char jobz = do_eigenvectors ? 'V' : 'N';
    const size_type n = A.nrows();
    const size_type lda = n; 

    // get optimal work size
    Vector<real_t<coeff_t>> w(n);
    size_type lwork = -1;
    std::vector<coeff_t> work(1);
    int info = 0;
    blaslapack::syev(&jobz, &uplo, &n, A.data(), &lda, w.data(), 
		     work.data(), &lwork, &info);
    assert(info == 0);
    lwork = static_cast<size_type>(real(work[0]));
    work.resize(lwork);
    
    // Run eigenvalue computation
    blaslapack::syev(&jobz, &uplo, &n, A.data(), &lda, w.data(), 
		     work.data(), &lwork, &info);
    assert(info == 0);
    return w;
  }

  template <class coeff_t>
  inline Vector<real_t<coeff_t>> 
  EigenvaluesH(const Matrix<coeff_t>& A, char uplo='U') {
    Matrix<coeff_t> A_copy = A;
    return EigenH(A_copy, false, uplo);
  }

}
#endif
