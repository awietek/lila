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

#ifndef LILA_SPARSE_LANCZOS_H_
#define LILA_SPARSE_LANCZOS_H_

#include <iostream>
#include <algorithm>
#include <cstdlib>

#include "../matrix.h"
#include "../tmatrix.h"
#include "../vector.h"
#include "../common.h"
#include "lanczos_convergence.h"

namespace lila {
    
  template <class vector_t>
  struct LanczosResults
  {
    using coeff_t = typename vector_t::coeff_type;
    using real_type = real_t<coeff_t>;
    
    Vector<real_type> eigenvalues;
    Tmatrix<real_type> tmatrix;
    real_type beta;
    
    std::vector<vector_t> vectors;
  };

  template <class vector_t, class multiply_f, class convergence_f>
  LanczosResults<vector_t>
  Lanczos(multiply_f A, vector_t& v0, convergence_f converged, 
	  const std::vector<Vector<typename vector_t::coeff_type>>& 
	  linear_combinations = {})
  {
    using coeff_type = typename vector_t::coeff_type;
    using real_type = real_t<coeff_type>;

    // Initialize Lanczos vectors and tmatrix
    vector_t v1 = v0;
    vector_t w = v0;

    Zeros(v0);
    Zeros(w);
    real_type alpha = 0; 
    real_type beta = 0;
    LanczosResults<vector_t> res = 
      { Vector<real_type>(), Tmatrix<real_type>(), 0, std::vector<vector_t>()};

    // Initialize linear combination vectors (e.g. eigenvectors / time evo)
    for (int i=0; i<(int)linear_combinations.size(); ++i)
      res.vectors.push_back(v0);

    // Normalize start vector or return if norm is zero
    double v1_norm = Norm(v1);
    if (!close(v1_norm, (real_type)0.)) 
      Normalize(v1);
    else
      {
	res.tmatrix = Tmatrix<real_type>(0);
	res.eigenvalues = Eigenvalues(res.tmatrix);
	res.beta = 0;
	return res;
      }

    // Main Lanczos loop
    int iteration = 0;
    while (!converged(res.tmatrix, beta))
      {
	// Build linear combinations
	int lin_combo_idx = 0;
	for (auto lin_combo : linear_combinations)
	  {
	    assert(iteration < lin_combo.size());
	    Add(v1, res.vectors[lin_combo_idx], lin_combo(iteration));
	    ++lin_combo_idx;
	  }
  
	// Lanczos recursion
	A(v1, w);                  // MVM
	alpha = real(Dot(v1, w));
	Add(v1, w, -(coeff_type)alpha);      // w -= alpha*v1;
	Add(v0, w, -(coeff_type)beta);       // w -= beta*v0;
	v0 = v1;
	v1 = w;

	// Extend Tmatrix
	if (iteration == 0) 
	  res.tmatrix = Tmatrix<real_type>(alpha);
	else 
	  res.tmatrix.push_back(alpha, beta);

	beta = Norm(v1);

	// Finish if Lanczos sequence is exhausted
	if (!close(beta, (real_type)0.)) 
	  Normalize(v1);
	else break;

	++iteration;
      }

    res.eigenvalues = Eigenvalues(res.tmatrix);
    res.beta = beta;
    return res;
  }

  template <class vector_t, class multiply_f>
  LanczosResults<vector_t>
  LanczosEigenvalues(multiply_f A, vector_t& v0,
		     real_t<typename vector_t::coeff_type> precision = 1e-12,
		     int n_eigenvalue = 0, 
		     std::string criterion = "Eigenvalues")
  {
    using real_type = real_t<typename vector_t::coeff_type>;
    using tmatrix_t = Tmatrix<real_type>;

    LanczosResults<vector_t> res;

    if (criterion == "Eigenvalues")
      {
	auto converged = 
	  [n_eigenvalue, precision](const tmatrix_t& tmat, real_type beta) { 
	    return LanczosConvergedEigenvalues(tmat, beta, n_eigenvalue, 
					       precision);
	  };
	res = Lanczos(A, v0, converged);
      }

    else if (criterion == "Ritz")
      {
	auto converged = 
	  [n_eigenvalue, precision](const tmatrix_t& tmat, real_type beta) {
	    return LanczosConvergedRitz(tmat, beta, n_eigenvalue, 
					precision);
	};
       res = Lanczos(A, v0, converged);
      }

    else
      {
	std::cerr << "Unknown convergence criterion for Lanczos algorithm!" 
		  << std::endl;
	exit(EXIT_FAILURE);
      }
    return res;
  }

  template <class vector_t, class multiply_f, class gen_t>
  LanczosResults<vector_t>
  LanczosEigenvectors(multiply_f A, vector_t v0, 
		      gen_t& gen, bool alter_generator,
		      real_t<typename vector_t::coeff_type> precision = 1e-12,
		      std::vector<int> num_eigenvectors = {0}, 
		      std::string criterion = "Ritz")
  {
    using coeff_type = typename vector_t::coeff_type;
    using real_type = real_t<coeff_type>;
    using tmatrix_t = Tmatrix<real_type>;

    // First run for Tmatrix
    Random(v0, gen, false);  // do not alter the generator
    int n_eigenvalue = 0;
    for (int i=0; i<(int)num_eigenvectors.size(); ++i)
      {
	assert(num_eigenvectors[i] >= 0);
	n_eigenvalue = std::max(n_eigenvalue, num_eigenvectors[i]);
      }

    auto res_first = 
      LanczosEigenvalues(A, v0, precision, n_eigenvalue, criterion);
    
    // Define fixed iterations convergence criterion
    int n_iterations = res_first.tmatrix.size();
    auto converged = 
      [n_iterations](const tmatrix_t& tmat, real_type beta) {
      (void)beta;
      return LanczosConvergedFixed(tmat, n_iterations);
    };
    
    // Prepare linear combinations
    auto evecs = Eigen(res_first.tmatrix).eigenvectors;
    std::vector<Vector<coeff_type>> linear_combinations;
    for (auto num_eigenvector : num_eigenvectors)
      {
	Vector<coeff_type> evc(n_iterations);
	for (int i=0; i<n_iterations; ++i)
	  evc(i) = (coeff_type)evecs(i, num_eigenvector);
	linear_combinations.push_back(evc);
      }

    // reset initial vector and rerun with linear combination
    Random(v0, gen, alter_generator);
    return Lanczos(A, v0, converged, linear_combinations);

  }
}
#endif
