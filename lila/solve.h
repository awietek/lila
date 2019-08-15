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

#ifndef LILA_SOLVE_H_
#define LILA_SOLVE_H_

#include "matrix.h"
#include "vector.h"
#include "blaslapack/blaslapack.h"
#include "special.h" 

namespace lila {

  /*! @brief Solves a system of linear equations with matrix r.h.s

    Performs the following operation:
    \f$ A^{-1}B \rightarrow B \f$

    @param A lila::Matrix (overwritten, contains factors L, U LU factorization) 
    @param B lila::Matrix (overwritten, contains solution matrix A^{-1}B)

    @tparam coeff_t type of coefficients of object 
   */
  template <class coeff_t>
  inline std::vector<int> Solve(Matrix<coeff_t>& A, Matrix<coeff_t>& B)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == A.ncols());
    assert(A.nrows() == B.nrows());
    const size_type n = A.nrows();
    const size_type n_rhs = B.ncols();
    const size_type lda = n;
    const size_type ldb = n;

    std::vector<int> ipiv(n); 
    int info = 0;
    blaslapack::gesv(&n, &n_rhs, A.data(), &lda, ipiv.data(), B.data(), &ldb, 
		     &info);
    return ipiv;
  }

  /*! @brief Solves a system of linear equations with vector r.h.s

    Performs the following operation:
    \f$ A^{-1}b \rightarrow b \f$

    @param A lila::Matrix (overwritten, contains factors L, U LU factorization) 
    @param B lila::Vector (overwritten, contains solution matrix A^{-1}b)

    @tparam coeff_t type of coefficients of object 
   */
  template <class coeff_t>
  inline std::vector<int> Solve(Matrix<coeff_t>& A, Vector<coeff_t>& X)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == X.nrows());
    const size_type n = A.nrows();
    const size_type n_rhs = 1;
    const size_type lda = n;
    const size_type ldb = n;

    std::vector<int> ipiv(n); 
    int info = 0;
    blaslapack::gesv(&n, &n_rhs, A.data(), &lda, ipiv.data(), X.data(), &ldb, 
		     &info);
    assert(info == 0);

    return ipiv;
  }

  template <class coeff_t>
  inline std::vector<int> LUDecompose(Matrix<coeff_t>& A)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    const size_type m = A.nrows();
    const size_type n = A.ncols();
    const size_type lda = m;

    std::vector<int> ipiv(n); 
    int info = 0;
    blaslapack::getrf(&m, &n, A.data(), &lda, ipiv.data(), &info);
    assert(info == 0);

    return ipiv;
  }

  template <class coeff_t>
  inline void LUSolve(const Matrix<coeff_t>& A, const std::vector<int>& ipiv,
		      Matrix<coeff_t>& B, char trans='N')
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == A.ncols());
    assert(A.nrows() == B.nrows());
    const size_type n = A.nrows();
    const size_type n_rhs = B.ncols();
    const size_type lda = n;
    const size_type ldb = n;

    int info = 0;
    blaslapack::getrs(&trans, &n, &n_rhs, A.data(), &lda, ipiv.data(), B.data(),
		      &ldb, &info);
    assert(info == 0);
  }

  template <class coeff_t>
  inline void LUSolve(const Matrix<coeff_t>& A, const std::vector<int>& ipiv,
		      Vector<coeff_t>& X, char trans='N')
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == X.nrows());
    const size_type n = A.nrows();
    const size_type n_rhs = 1;
    const size_type lda = n;
    const size_type ldb = n;

    int info = 0;
    blaslapack::getrs(&trans, &n, &n_rhs, A.data(), &lda, ipiv.data(), X.data(),
		      &ldb, &info);
    assert(info == 0);
  }
  
  template <class coeff_t>
  inline void Invert(Matrix<coeff_t>& A)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == A.ncols());
    const size_type n = A.nrows();

    std::vector<int> ipiv(n); 
    size_type lwork = n*n;
    std::vector<coeff_t> work(lwork);
    int info = 0;

    blaslapack::getrf(&n, &n, A.data(), &n, ipiv.data(), &info);
    assert(info == 0);
    blaslapack::getri(&n, A.data(), &n, ipiv.data(), work.data(), &lwork, 
		      &info);
    assert(info == 0);
  }

  template <class coeff_t>
  inline std::vector<coeff_t> QRDecompose(Matrix<coeff_t>& A)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    const size_type m = A.nrows();
    const size_type n = A.ncols();
    const size_type lda = m;

    std::vector<coeff_t> tau(std::min(m, n)); 
    int info = 0;
    
    // get optimal work size
    size_type lwork = -1;
    std::vector<coeff_t> work(1);
    blaslapack::geqrf(&m, &n, A.data(), &lda, tau.data(), work.data(), &lwork, 
		      &info);
    assert(info == 0);
    lwork = static_cast<size_type>(real(work[0]));
    work.resize(lwork);

    // Run QR Decomposition
    blaslapack::geqrf(&m, &n, A.data(), &lda, tau.data(), work.data(), &lwork, 
		      &info);
    assert(info == 0);

    return tau;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> QRGetQ(Matrix<coeff_t>& A, 
				const std::vector<coeff_t>& tau)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    const size_type m = A.nrows();
    const size_type k = tau.size();
    const size_type lda = m;

    Matrix<coeff_t> Q(A);
    int info = 0;

    // get optimal work size
    int lwork = -1;
    std::vector<coeff_t> work(1);
    blaslapack::orgqr(&m, &k, &k, Q.data(), &lda, tau.data(), work.data(), 
		      &lwork, &info);
    assert(info == 0);
    lwork = static_cast<size_type>(real(work[0]));
    work.resize(lwork);

    // Compute Q
    blaslapack::orgqr(&m, &k, &k, Q.data(), &lda, tau.data(), work.data(), 
		      &lwork, &info);

    Q.resize(m, k);
    assert(info == 0);
    return Q;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> GetUpper(Matrix<coeff_t>& A)
  {
    const size_type m = A.nrows();
    const size_type n = A.ncols();
    const size_type k = std::min(m, n);
    
    Matrix<coeff_t> R(k, n);
    Zeros(R);
    for (int row=0; row < k; ++row)
      for (int column=0; column < n; ++column)
	if (row <= column) R(row, column) = A(row, column);
    return R;			     
  }

}

#endif
