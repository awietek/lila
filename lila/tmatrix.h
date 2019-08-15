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

#ifndef LILA_SPARSE_TMATRIX_H_
#define LILA_SPARSE_TMATRIX_H_

#include "vector.h"
#include "matrix.h"
#include "special.h"
#include "eigen.h"

namespace lila {


  /*! @brief Class describing real, tridiagonal, symmetric matrices

    __Usage example __
    @code
    int seed = 42;
    uniform_dist_t<coeff_t> fdist(-1., 1.);
    uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    Vector<coeff_t> diag(n);
    Vector<coeff_t> offdiag(n-1);
    Random(diag, fgen);
    Random(offdiag, fgen);

    Tmatrix<coeff_t> tmat(diag, offdiag);
    auto fullmat = Matrix<coeff_t>(tmat);
    auto eigs = Eigenvalues(tmat);
    auto res = Eigen(tmat);

    LilaPrint(tmat.diag());
    LilaPrint(tmat.offdiag());
    LilaPrint(fullmat);
    LilaPrint(eigs);

    auto evecs = res.eigenvectors;
    for (int i : evecs.cols())
    {
      printf("norm: %f\n", Dot(evecs.col(i), evecs.col(i)));
      printf("eig : %f\n", eigs(i));
      coeff_t test = Dot(evecs.col(i), Mult(fullmat, evecs.col(i)));
      printf("test: %f\n", test);
    }
    @endcode
   */
  template <class coeff_t>
  class Tmatrix 
  {
  public:
    using size_type = lila::size_type; 
    using coeff_type = coeff_t;
    using value_type = coeff_t;        
    using vector_type = std::vector<coeff_t>;

    Tmatrix() = default;
    Tmatrix(const Vector<coeff_t>& diag, 
	    const Vector<coeff_t>& offdiag)
      : diag_(diag), 
	offdiag_(offdiag)
    { assert(diag.size() == offdiag.size() + 1); }

    Tmatrix(const coeff_t& alpha0)
    { diag_.push_back(alpha0); }
    
    size_type size() const { return diag_.size(); }
    Vector<coeff_t> diag() const { return diag_; }
    Vector<coeff_t> offdiag() const { return offdiag_; }
    void push_back(const coeff_t& alpha, const coeff_t& beta)
    {
      diag_.push_back(alpha);
      offdiag_.push_back(beta);
    }

    void resize(size_type size) {
      diag_.resize(size);
      if (size > 0)
	offdiag_.resize(size-1);
      else
	offdiag_.resize(0);
    }

    /// conversion to a regular lila::Matrix
    operator Matrix<coeff_t>() const
    { 
      return Diag<coeff_t>(diag_) + Diag<coeff_t>(offdiag_, -1) + 
	Diag<coeff_t>(offdiag_, 1);
    }
    

  private:
    Vector<coeff_t> diag_;
    Vector<coeff_t> offdiag_;
  };
  
  /// Calculates eigenvalues of T matrix
  template <class coeff_t>
  inline Vector<coeff_t> Eigenvalues(const Tmatrix<coeff_t>& tmat)
  {
    Vector<coeff_t> eigs = tmat.diag();
    Vector<coeff_t> offdiag = tmat.offdiag();

    assert(eigs.size() == offdiag.size() + 1);
      
    char job = 'N';
    int N = eigs.size();
    int ldz = N;
    int info = 0;
    blaslapack::stev(&job, &N, eigs.data(), offdiag.data(), 
		     NULL, &ldz, NULL, &info);
    return eigs;
  }


  /// Calculates eigenvalues and eigenvectors of T matrix
  template <class coeff_t>
  inline EigenSymResults<coeff_t> Eigen(const Tmatrix<coeff_t>& tmat)
  {
    EigenSymResults<coeff_t> results;
    results.eigenvalues = tmat.diag();
    Vector<coeff_t> offdiag = tmat.offdiag();
    
    assert(results.eigenvalues.size() == offdiag.size() + 1);

    char job = 'V';
    int N = results.eigenvalues.size();
    int ldz = N;
    int info = 0;

    results.eigenvectors = Matrix<coeff_t>(N, N);
    
    Vector<coeff_t> work(2*N - 2);
    blaslapack::stev(&job, &N, results.eigenvalues.data(), offdiag.data(),
		     results.eigenvectors.data(), &ldz, work.data(), &info);
    return results;
  }
    

  template <class coeff_t>
  inline void PrintPretty(const char* identifier, 
  			  const Tmatrix<coeff_t>& tmat) {
    printf("%s:\n", identifier);
    PrintPretty("diag", tmat.diag());
    PrintPretty("offdiag", tmat.offdiag());
  }

}

#endif
