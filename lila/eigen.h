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
#include "print.h"

namespace lila {


  /*! @brief Struct holding results of symmetric (or hermitian) eigenvalue
             problems
   */
  template <class coeff_t>
  struct EigenSymResults {
    Vector<real_t<coeff_t>> eigenvalues;
    Matrix<coeff_t> eigenvectors;
  };

  template <class coeff_t>
  inline Vector<real_t<coeff_t>> 
  EigenSymInplace(Matrix<coeff_t>& A, bool do_eigenvectors=true, char uplo='U')
  {
    using size_type = blaslapack::blas_size_t;
    assert(A.nrows() == A.ncols());

    char jobz = do_eigenvectors ? 'V' : 'N';
    size_type n = A.nrows();
    size_type lda = n; 

    // get optimal work size
    Vector<real_t<coeff_t>> w(n);
    size_type lwork = -1;
    std::vector<coeff_t> work(1);
    int info = 0;
    
    blaslapack::syev(&jobz,
		     &uplo,
		     &n,
		     LILA_BLAS_CAST(coeff_t,A.data()),
		     &lda,
		     LILA_BLAS_CAST(real_t<coeff_t>,w.data()), 
		     LILA_BLAS_CAST(coeff_t,work.data()),
		     &lwork,
		     &info);
    
    assert(info == 0);
    lwork = static_cast<size_type>(real(work[0]));
    work.resize(lwork);
    
    // Run eigenvalue computation
    blaslapack::syev(&jobz,
		     &uplo,
		     &n,
		     LILA_BLAS_CAST(coeff_t,A.data()),
		     &lda,
		     LILA_BLAS_CAST(real_t<coeff_t>,w.data()), 
		     LILA_BLAS_CAST(coeff_t,work.data()),
		     &lwork,
		     &info);
    
    assert(info == 0);
    return w;
  }

  template <class coeff_t>
  inline EigenSymResults<coeff_t>
  EigenSym(const Matrix<coeff_t>& A, char uplo='U')
  {
    EigenSymResults<coeff_t> result;
    result.eigenvectors = A;
    result.eigenvalues = EigenSymInplace(result.eigenvectors, true, uplo);
    return result;  
  }

  template <class coeff_t>
  inline Vector<real_t<coeff_t>> 
  EigenvaluesSym(const Matrix<coeff_t>& A, char uplo='U') {
    Matrix<coeff_t> A_copy = A;
    return EigenSymInplace(A_copy, false, uplo);
  }



  template <class coeff_t>
  inline Vector<real_t<coeff_t>> 
  EigenGenSymDefInplace(Matrix<coeff_t>& A, Matrix<coeff_t>& B, 
			bool do_eigenvectors=true, char uplo='U', 
			int itype = 1)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == A.ncols());
    assert(B.nrows() == B.ncols());
    assert(A.nrows() == B.nrows());

    char jobz = do_eigenvectors ? 'V' : 'N';
    size_type n = A.nrows();
    size_type lda = n;
    size_type ldb = n;
    int info = 0;

    // get optimal work size
    Vector<real_t<coeff_t>> w(n);
    size_type lwork = -1;
    std::vector<coeff_t> work(1);
    blaslapack::syev(&jobz,
		     &uplo,
		     &n,
		     LILA_BLAS_CAST(coeff_t,A.data()),
		     &lda,
		     LILA_BLAS_CAST(real_t<coeff_t>,w.data()), 
		     LILA_BLAS_CAST(coeff_t,work.data()),
		     &lwork,
		     &info);
    assert(info == 0);
    lwork = static_cast<size_type>(real(work[0]));
    work.resize(lwork);

    // Run eigenvalue computation
    blaslapack::sygv(&itype,
		     &jobz,
		     &uplo,
		     &n,
		     LILA_BLAS_CAST(coeff_t,A.data()),
		     &lda,
		     LILA_BLAS_CAST(coeff_t,B.data()),
		     &ldb,
		     LILA_BLAS_CAST(real_t<coeff_t>,w.data()),
		     LILA_BLAS_CAST(coeff_t,work.data()),
		     &lwork,
		     &info);

    assert(info == 0);
    return w;
  }

  template <class coeff_t>
  inline EigenSymResults<coeff_t>
  EigenGenSymDef(const Matrix<coeff_t>& A, const Matrix<coeff_t>& B,
		 char uplo='U', int itype = 1)
  {
    EigenSymResults<coeff_t> result;
    result.eigenvectors = A;
    Matrix<coeff_t> B_copy = B;
    result.eigenvalues = EigenGenSymDefInplace(result.eigenvectors, B_copy,
					       true, uplo, itype);
    return result;  
  }

  template <class coeff_t>
  inline Vector<real_t<coeff_t>>
  EigenvaluesGenSymDef(const Matrix<coeff_t>& A, const Matrix<coeff_t>& B,
		     char uplo='U', int itype = 1)
  {
    Matrix<coeff_t> A_copy = A;
    Matrix<coeff_t> B_copy = B;
    return EigenGenSymDefInplace(A_copy, B_copy, false, uplo, itype);

  }



  
  template <class coeff_t>
  struct EigenResults {
    Vector<complex_t<coeff_t>> eigenvalues;
    Matrix<coeff_t> left_eigenvectors;
    Matrix<coeff_t> right_eigenvectors;    
  };

  template <class coeff_t>
  inline EigenResults<coeff_t> 
  Eigen(Matrix<coeff_t>& A, bool do_right_eigenvectors=true, 
	bool do_left_eigenvectors=true)
  {
    using size_type = blaslapack::blas_size_t;
    assert(A.nrows() == A.ncols());

    char jobvl = do_left_eigenvectors ? 'V' : 'N';
    char jobvr = do_right_eigenvectors ? 'V' : 'N';
    size_type n = A.nrows();
    size_type lda = n;
 
    size_type ldvl = do_left_eigenvectors ? n : 1;
    size_type ldvr = do_right_eigenvectors ? n : 1;
    
    EigenResults<coeff_t> result;
    result.eigenvalues.resize(n);
    if (do_left_eigenvectors) result.left_eigenvectors.resize(ldvl, n);
    if (do_right_eigenvectors) result.right_eigenvectors.resize(ldvr, n);

    // get optimal work size
    Vector<complex_t<coeff_t>>& w = result.eigenvalues;
    size_type lwork = -1;
    std::vector<coeff_t> work(1);
    int info = 0;
    blaslapack::geev(&jobvl,
		     &jobvr,
		     &n,
		     LILA_BLAS_CAST(coeff_t,A.data()),
		     &lda,
		     LILA_BLAS_CAST(complex_t<coeff_t>,w.data()), 
		     LILA_BLAS_CAST(coeff_t, result.left_eigenvectors.data()),
		     &ldvl,
		     LILA_BLAS_CAST(coeff_t, result.right_eigenvectors.data()),
		     &ldvr, 
		     LILA_BLAS_CAST(coeff_t, work.data()),
		     &lwork,
		     &info);
    assert(info == 0);
    lwork = static_cast<size_type>(real(work[0]));
    work.resize(lwork);

    // Run eigenvalue computation
    blaslapack::geev(&jobvl,
		     &jobvr,
		     &n,
		     LILA_BLAS_CAST(coeff_t,A.data()),
		     &lda,
		     LILA_BLAS_CAST(complex_t<coeff_t>,w.data()), 
		     LILA_BLAS_CAST(coeff_t,result.left_eigenvectors.data()),
		     &ldvl,
		     LILA_BLAS_CAST(coeff_t,result.right_eigenvectors.data()),
		     &ldvr, 
		     LILA_BLAS_CAST(coeff_t,work.data()),
		     &lwork,
		     &info);
    assert(info == 0);
    return result;
  }

  template <class coeff_t>
  inline Vector<complex_t<coeff_t>> 
  Eigenvalues(Matrix<coeff_t>& A) {
    EigenResults<coeff_t> res = Eigen(A);
    return res.eigenvalues;
  }

}
#endif
