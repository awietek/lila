// Copyright 2019 Alexander Wietek - All Rights Reserved.
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

#ifndef LILA_SPARSE_LANCZOSCONVERGENCE_H_
#define LILA_SPARSE_LANCZOSCONVERGENCE_H_

#include "../tmatrix.h"
#include "../vector.h"
#include "../common.h"

namespace lila {

  template <class coeff_t>
  bool LanczosConvergedEigenvalues
  (const Tmatrix<coeff_t>& tmat, coeff_t beta, int n_eigenvalue, 
   coeff_t precision)
  {
    int size = tmat.size();
    if (size < 2) return false;
    else
      {
	if (close(beta, (coeff_t)0.)) return true;
	if (size <= n_eigenvalue) return false;

	// Compute previous tmatrix
	auto tmat_previous = tmat;
	tmat_previous.resize(size - 1);

	auto eigs = Eigenvalues(tmat);
	auto eigs_previous = Eigenvalues(tmat_previous);

	double residue = 
	  std::abs(eigs(n_eigenvalue) - eigs_previous(n_eigenvalue)) /
	  std::abs(eigs(n_eigenvalue));

	return (residue < precision);
      }
  }
  

  template <class coeff_t>
  bool LanczosConvergedRitz
  (const Tmatrix<coeff_t>& tmat, coeff_t beta, int n_eigenvalue, 
   coeff_t precision)
  {
    int size = tmat.size();
    if (size < 2) return false;
    else
      {
	// Lanczos sequence exhausted
	if (close(beta, (coeff_t)0.)) return true;
	if (size <= n_eigenvalue) return false;

	// Compute all Ritz residuals
	auto evecs = Eigen(tmat).eigenvectors;
	bool conv = true;
	for (int n = 0; n <= n_eigenvalue; ++n)
	  {
	    int idx = size - 1 - n_eigenvalue;
	    double residue = std::abs(evecs(idx, size-1) * beta);
	    conv &= (residue < precision);
	  }
	return conv;
      }
  }


  template <class coeff_t>
  bool LanczosConvergedFixed
  (const Tmatrix<coeff_t>& tmat, int n_iterations)
  {
    int size = tmat.size();
    return size >= n_iterations;
  }

  template <class coeff_t, class ccoeff_t>
  bool LanczosConvergedTimeEvolution
  (const Tmatrix<coeff_t>& tmat, coeff_t beta, ccoeff_t tau, 
   coeff_t precision, int max_iter, coeff_t nrm)
  {
    int size = tmat.size();
    if (size < 2) return false;
    else
      {
	// Lanczos sequence exhausted
	if (close(beta, (coeff_t)0.)) return true;

	// Prepare extended T-matrix for exponentiation
	auto tmatr = Matrix<coeff_t>(tmat);
	auto tmat_ext = Zeros<ccoeff_t>(size, size);
	for (auto i : tmatr.rows())
	  for (auto j : tmatr.cols())
	    tmat_ext(i,j) = (ccoeff_t)tmatr(i,j);

	int ext_size = size + 2;
	tmat_ext.resize(ext_size, ext_size);
	
	tmat_ext(size-1, size) = beta;
	tmat_ext(size, size-1) = beta;
	
	tmat_ext(size-1, size) = 0.;
	tmat_ext(size+1, size) = 1.;
    
	// Exponentiate extended T-matrix
	auto tmat_ext_exp = ExpM(tmat_ext, tau);
	coeff_t phi1 = std::abs( nrm * tmat_ext_exp(size, 0) );
	coeff_t phi2 = std::abs( nrm * tmat_ext_exp(size + 1, 0));
	// printf("phi1: %g, phi2: %g\n", phi1, phi2);
	coeff_t error;      
	if (phi1 > 10*phi2) error = phi2;
	else if (phi1 > phi2) error = (phi1*phi2)/(phi1-phi2);
	else error = phi1;
	if ((error < precision) || (size == max_iter-1))
	  {
	    if (size == max_iter-1) 
	      printf("warning: lanczosTevol not converged in %d steps\n", max_iter);
	    return true;
	  }
	else return false;
      }
  }
      
}  // namespace lila

#endif
