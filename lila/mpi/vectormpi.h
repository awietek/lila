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

#ifndef LILA_MPI_VECTORMPI_H_
#define LILA_MPI_VECTORMPI_H_

#include <mpi.h>

#include "../vector.h"
#include "../common.h"
#include "stable_dot_product.h"

namespace lila {
  
  template <class coeff_t>
  class VectorMPI {
  public:
    using size_type = int64; 
    using coeff_type = coeff_t; 
    using value_type = coeff_t;   

    VectorMPI() : size_(0), vector_local_(0)
    { }
    ~VectorMPI() = default;
    VectorMPI(const VectorMPI&) = default;
    VectorMPI& operator=(VectorMPI&) = default; 
    VectorMPI(VectorMPI&&) = default;
    VectorMPI& operator=(VectorMPI&&) = default;

    explicit VectorMPI(size_type size_local) 
      : size_(0), vector_local_(size_local) 
    { MPI_Allreduce(&size_local, &size_, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD); }
    
    explicit VectorMPI(Vector<coeff_t>& vec) 
      : size_(0), vector_local_(vec) 
    { 
      int64 size_local = (int64)vec.size();
      MPI_Allreduce(&size_local, &size_, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD); 
    }

    coeff_t operator()(size_type i) const { return vector_local_(i); }
    coeff_t& operator()(size_type i) { return vector_local_(i); }

    void clear() { 
      size_ = 0;
      vector_local_.clear(); 
    }
    void shrink_to_fit() { vector_local_.shrink_to_fit(); }

    size_type size() const { return vector_local_.size(); }
    size_type size_global() const { return size_; }
    void resize(size_type size_local)
    {
      vector_local_.resize(size_local);
      MPI_Allreduce(&size_local, &size_, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);  
    }

    lila::Vector<coeff_t>&  vector_local() { return vector_local_; }
    const lila::Vector<coeff_t>&  vector_local() const 
    { return vector_local_; }

    coeff_t* data() { return vector_local_.data(); }
    const coeff_t* data() const { return vector_local_.data(); }
    
  private:
    uint64 size_;
    lila::Vector<coeff_t> vector_local_;
  };


  template <class coeff_t, class function_t>
  inline VectorMPI<coeff_t> Map(const VectorMPI<coeff_t>& X, function_t func)
  { 
    std::for_each(X.data(), X.data() + X.size(), func);
    return X;
  }


  template <class coeff_t>
  inline void Copy(const VectorMPI<coeff_t>& X,  VectorMPI<coeff_t>& Y)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = (size_type)X.size();
    const size_type dy = (size_type)Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::copy(&dx, X.data(), &inc, Y.data(), &inc);
  }


  template <class coeff_t>
  inline void Add(const VectorMPI<coeff_t>& X,  VectorMPI<coeff_t>& Y, 
		  coeff_t alpha = static_cast<coeff_t>(1.))
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = (size_type)X.size();
    const size_type dy = (size_type)Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    blaslapack::axpy(&dx, &alpha, X.data(), &inc, Y.data(), &inc);
  }


  template <class coeff_t>
  inline void Scale(const coeff_t& alpha,  VectorMPI<coeff_t>& X)
  {
    using size_type = blaslapack::blas_size_t;
    const size_type dx = (size_type)X.size();
    const size_type inc = 1;
    blaslapack::scal(&dx, &alpha, X.data(), &inc);
  }


  template <class coeff_t>
  inline coeff_t Dot(const VectorMPI<coeff_t>& X, const VectorMPI<coeff_t>& Y)
  {
    const uint64 dx = (uint64)X.size();
    const uint64 dy = (uint64)Y.size();
    assert(dx == dy); // Check if valid dimensions
    const size_type inc = 1;
    return stable_dot_product(dx, X.data(), Y.data());
  }


  template <class coeff_t>
  inline real_t<coeff_t> Norm(const VectorMPI<coeff_t>& X)
  { return sqrt(real(Dot(X, X))); }

  template <class coeff_t, class gen_t>
  void Random(VectorMPI<coeff_t>& vec, gen_t& gen, bool alter_generator=true)
  { 
    if(alter_generator) std::for_each(vec.vector_local().begin(), 
				      vec.vector_local().end(), 
				      [&gen](coeff_t& c){ c = gen(); });
    else std::generate(vec.vector_local().begin(), vec.vector_local().end(), gen); 
  }

}

#endif
