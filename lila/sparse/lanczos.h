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

#include "../matrix.h"
#include "tmatrix.h"
#include "../vector.h"
#include "../common.h"

namespace lila {
    
  template <class coeff_t, class multiply_f, class vector_t = Vector<coeff_t>>
  class Lanczos
  {

  public:
    Lanczos(uint64 dimension, int random_seed, int max_iterations, 
	    double precision, int num_eigenvalue, multiply_f multiply) 
      : dimension_(dimension),
	random_seed_(random_seed),
	max_iterations_(max_iterations),
	precision_(precision),
	num_eigenvalue_(num_eigenvalue),
	multiply_(multiply),
	init_state_(nullptr)
    { }

    Vector<double> alphas() const { return alphas_; }
    Vector<double> betas() const { return betas_; }
    Tmatrix<double> tmatrix() const { 
      // remove last element from betas_
      Vector<double> betas_short = betas_;
      betas_short.resize(betas_.size()-1);
      return Tmatrix<double>(alphas_, betas_short); 
    }

    void set_init_state(vector_t& init_state)
    { 
      assert((uint64)init_state.size() == dimension_);
      init_state_ = &init_state; 
    }

    bool finished() const 
    { 
      if (alphas_.size() > max_iterations_) return true;
      else return converged();
    }

    bool converged() const 
    { 
      int iter = alphas_.size();
      if (iter <= std::max(2, num_eigenvalue_)) return false;
      else
	{
	  // // total norm distance
	  // Matrix<double> evecs = Eigen(tmatrix()).eigenvectors;
	  // double e = evecs(iter-1, num_eigenvalue_);
	  // double b = betas_(iter-1);
	  // if (close(b, 0.)) return true;
	  // double residue = std::abs(e * b); 


	  // eigenvalue convergence
	  auto tmat = tmatrix();
	  Vector<double> alphas_short = alphas_;
	  alphas_short.resize(alphas_.size()-1);
	  auto eigs = Eigenvalues(tmat);

	  Vector<double> betas_short = betas_;
	  betas_short.resize(betas_.size()-2);
	  auto prev_tmat = Tmatrix<double>(alphas_short, betas_short); 
	  auto prev_eigs = Eigenvalues(prev_tmat);
	  
	  double residue = std::abs(prev_eigs(num_eigenvalue_) - 
	  			    eigs(num_eigenvalue_));

	  return (residue < precision_);
	}
    }

    Vector<double> eigenvalues()
    {
      // Initialize Lanczos vectors
      vector_t v0(dimension_);
      vector_t v1(dimension_);
      vector_t w(dimension_);
      if (!init_state_)
	{
	  lila::normal_dist_t<real_t<coeff_t>> ddist(0., 1.);
	  lila::normal_gen_t<coeff_t> gen(ddist, random_seed_);
	  Random(v1, gen);
	}
      else v1 = *init_state_;
      // v1 /= (coeff_t)Norm(v1);
      coeff_t norm = Norm(v1);
      Scale((coeff_t)1./norm, v1);

      double alpha, beta;
      int iter = 0;

      // iterate and build T-matrix
      while (!finished())
	{	  
	  lanczos_step(v0, v1, w, alpha, beta, iter);
	  alphas_.push_back(alpha);
	  betas_.push_back(beta);
	  ++iter;
	}
      return Eigenvalues(tmatrix());
    }

    std::vector<vector_t> eigenvectors(std::vector<int> num_eigenvectors)
    {
      int tmat_dim = alphas_.size();
      assert(tmat_dim > 0);

      // Prepare eigenvector computation
      auto t_res = Eigen(tmatrix());
      auto t_evecs = t_res.eigenvectors;
      std::vector<vector_t> evecs;
      for (int k=0; k < (int)num_eigenvectors.size(); ++k)
	{
	  assert(num_eigenvectors[k] < tmat_dim);
	  evecs.push_back(vector_t(dimension_));
	}

      // Initialize Lanczos vectors
      vector_t v0(dimension_);
      vector_t v1(dimension_);
      vector_t w(dimension_);

      if (!init_state_)
	{
	  lila::normal_dist_t<real_t<coeff_t>> ddist(0., 1.);
	  lila::normal_gen_t<coeff_t> gen(ddist, random_seed_);
	  Random(v1, gen);
	}
      else v1 = *init_state_;
      // v1 /= (coeff_t)Norm(v1);
      coeff_t norm = Norm(v1);
      Scale((coeff_t)1./norm, v1);

      // Reiterate to construct eigenvectors
      double alpha, beta;
      for (int iter = 0; iter < tmat_dim; ++iter )
	{
	  // Build linear combinations
	  for (int k=0; k < (int)num_eigenvectors.size(); ++k)
	    {
	      double coeff = t_evecs(iter, num_eigenvectors[k]);
	      evecs[k] += (coeff_t)coeff * v1;
	    }	  
	  lanczos_step(v0, v1, w, alpha, beta, iter);
	}
      return evecs;
    }

  private:
    void lanczos_step(vector_t& v0, vector_t& v1, vector_t& w, double& alpha, 
		      double& beta, int iter)
    {
      beta = (iter==0) ? 0. : betas_(iter-1);
      multiply_(v1, w);
      alpha = (double)real(Dot(v1, w));

      // w -= (coeff_t)alpha*v1;
      Add(v1, w, -(coeff_t)alpha);

      // w -= (coeff_t)beta*v0;
      Add(v0, w, -(coeff_t)beta);

      v0 = v1;
      beta = Norm(w);
      v1 = w;

      if (!close(beta, 0.))
	{
	  // v1 /= (coeff_t)beta;
	  Scale((coeff_t)(1./beta), v1);
	}
    }
    
    uint64 dimension_;
    int random_seed_;
    int max_iterations_;
    double precision_;
    int num_eigenvalue_;

    multiply_f multiply_;
    Vector<double> alphas_;
    Vector<double> betas_;

    vector_t* init_state_;
  };
      
}

#endif
