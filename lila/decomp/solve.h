#pragma once

#include <lila/matrix.h>
#include <lila/vector.h>
#include <lila/special/special.h>
#include <lila/blaslapack/blaslapack.h>

namespace lila {

  template <class coeff_t>
  inline std::vector<int> Solve(Matrix<coeff_t>& A, Matrix<coeff_t>& B)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == A.ncols());
    assert(A.nrows() == B.nrows());
    size_type n = A.nrows();
    size_type n_rhs = B.ncols();
    size_type lda = n;
    size_type ldb = n;

    std::vector<int> ipiv(n); 
    int info = 0;
    blaslapack::gesv(&n,
		     &n_rhs,
		     LILA_BLAS_CAST(coeff_t,A.data()),
		     &lda,
		     ipiv.data(),
		     LILA_BLAS_CAST(coeff_t,B.data()),
		     &ldb, 
		     &info);
    return ipiv;
  }

  template <class coeff_t>
  inline std::vector<int> Solve(Matrix<coeff_t>& A, Vector<coeff_t>& X)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == X.nrows());
    size_type n = A.nrows();
    size_type n_rhs = 1;
    size_type lda = n;
    size_type ldb = n;

    std::vector<int> ipiv(n); 
    int info = 0;
    blaslapack::gesv(&n,
		     &n_rhs,
		     LILA_BLAS_CAST(coeff_t,A.data()),
		     &lda,
		     ipiv.data(),
		     LILA_BLAS_CAST(coeff_t,X.data()),
		     &ldb, 
		     &info);
    assert(info == 0);

    return ipiv;
  }

  template <class coeff_t>
  inline std::vector<int> LUDecompose(Matrix<coeff_t>& A)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    size_type m = A.nrows();
    size_type n = A.ncols();
    size_type lda = m;

    std::vector<int> ipiv(n); 
    int info = 0;
    blaslapack::getrf(&m,
		      &n,
		      LILA_BLAS_CAST(coeff_t,A.data()),
		      &lda,
		      ipiv.data(),
		      &info);
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
    size_type n = A.nrows();
    size_type n_rhs = B.ncols();
    size_type lda = n;
    size_type ldb = n;

    int info = 0;
    blaslapack::getrs(&trans,
		      &n,
		      &n_rhs,
		      LILA_BLAS_CONST_CAST(coeff_t,A.data()),
		      &lda,
		      LILA_BLAS_CONST_CAST(int, ipiv.data()),
		      LILA_BLAS_CAST(coeff_t,B.data()),
		      &ldb,
		      &info);
    assert(info == 0);
  }

  template <class coeff_t>
  inline void LUSolve(const Matrix<coeff_t>& A, const std::vector<int>& ipiv,
		      Vector<coeff_t>& X, char trans='N')
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == X.nrows());
    size_type n = A.nrows();
    size_type n_rhs = 1;
    size_type lda = n;
    size_type ldb = n;

    int info = 0;
    blaslapack::getrs(&trans,
		      &n,
		      &n_rhs,
		      LILA_BLAS_CONST_CAST(coeff_t, A.data()),
		      &lda,
		      LILA_BLAS_CONST_CAST(int, ipiv.data()),
		      LILA_BLAS_CAST(coeff_t, X.data()),
		      &ldb,
		      &info);
    assert(info == 0);
  }
  
  template <class coeff_t>
  inline void Invert(Matrix<coeff_t>& A)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == A.ncols());
    size_type n = A.nrows();

    std::vector<int> ipiv(n); 
    size_type lwork = n*n;
    std::vector<coeff_t> work(lwork);
    int info = 0;

    blaslapack::getrf(&n,
		      &n,
		      LILA_BLAS_CAST(coeff_t, A.data()),
		      &n,
		      ipiv.data(),
		      &info);
    assert(info == 0);
    blaslapack::getri(&n,
		      LILA_BLAS_CAST(coeff_t ,A.data()),
		      &n,
		      ipiv.data(),
		      LILA_BLAS_CAST(coeff_t, work.data()),
		      &lwork, 
		      &info);
    assert(info == 0);
  }

  template <class coeff_t>
  inline std::vector<coeff_t> QRDecompose(Matrix<coeff_t>& A)
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    size_type m = A.nrows();
    size_type n = A.ncols();
    size_type lda = m;

    std::vector<coeff_t> tau(std::min(m, n)); 
    int info = 0;
    
    // get optimal work size
    size_type lwork = -1;
    std::vector<coeff_t> work(1);
    blaslapack::geqrf(&m,
		      &n,
		      LILA_BLAS_CAST(coeff_t,A.data()),
		      &lda,
		      LILA_BLAS_CAST(coeff_t,tau.data()),
		      LILA_BLAS_CAST(coeff_t,work.data()),
		      &lwork, 
		      &info);
    assert(info == 0);
    lwork = static_cast<size_type>(real(work[0]));
    work.resize(lwork);

    // Run QR Decomposition
    blaslapack::geqrf(&m,
		      &n,
		      LILA_BLAS_CAST(coeff_t,A.data()),
		      &lda,
		      LILA_BLAS_CAST(coeff_t,tau.data()),
		      LILA_BLAS_CAST(coeff_t,work.data()),
		      &lwork, 
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
    size_type m = A.nrows();
    size_type k = tau.size();
    size_type lda = m;

    Matrix<coeff_t> Q(A);
    int info = 0;

    // get optimal work size
    int lwork = -1;
    std::vector<coeff_t> work(1);
    blaslapack::orgqr(&m,
		      &k,
		      &k,
		      LILA_BLAS_CAST(coeff_t,Q.data()),
		      &lda,
		      LILA_BLAS_CONST_CAST(coeff_t,tau.data()),
		      LILA_BLAS_CAST(coeff_t,work.data()), 
		      &lwork,
		      &info);
    assert(info == 0);
    lwork = static_cast<size_type>(real(work[0]));
    work.resize(lwork);

    // Compute Q
    blaslapack::orgqr(&m,
		      &k,
		      &k,
		      LILA_BLAS_CAST(coeff_t,Q.data()),
		      &lda,
		      LILA_BLAS_CONST_CAST(coeff_t,tau.data()),
		      LILA_BLAS_CAST(coeff_t,work.data()), 
		      &lwork,
		      &info);

    Q.resize(m, k);
    assert(info == 0);
    return Q;
  }

  template <class coeff_t>
  inline Matrix<coeff_t> GetUpper(Matrix<coeff_t>& A)
  {
    size_type m = A.nrows();
    size_type n = A.ncols();
    size_type k = std::min(m, n);
    
    Matrix<coeff_t> R(k, n);
    Zeros(R);
    for (int row=0; row < k; ++row)
      for (int column=0; column < n; ++column)
	if (row <= column) R(row, column) = A(row, column);
    return R;			     
  }

}
