#pragma once

#include <lila/matrix.h>
#include <lila/blaslapack/blaslapack.h>

namespace lila {

  template <class coeff_t>
  inline void 
  CholeskyInplace(Matrix<coeff_t>& A, char uplo = 'U')
  {
    using size_type = blaslapack::blas_size_t;

    // check / get dimensions
    assert(A.nrows() == A.ncols());
    size_type n = A.nrows();
    size_type lda = n;
    int info = 0;
    blaslapack::potrf(&uplo,
		      &n,
		      LILA_BLAS_CAST(coeff_t,A.data()),
		      &lda,
		      &info);
  }

  template <class coeff_t>
  inline Matrix<coeff_t> 
  Cholesky(const Matrix<coeff_t>& A, char uplo = 'U')
  {
    auto res = A;
    CholeskyInplace(res, uplo);

    // Set lower triangular to zero
    int n = A.nrows();
    if (uplo == 'U')
      {
	for (int j=0; j<n; ++j)	
	  for (int i=j+1; i<n; ++i)
	    res(i, j) = 0;
      }
    else if (uplo == 'L')
      {
	for (int i=0; i<n; ++i)	
	  for (int j=i+1; j<n; ++j)
	    res(i, j) = 0;
      }
    return res;
  }


  

}
