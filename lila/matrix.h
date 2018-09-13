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

#ifndef LILA_MATRIX_H_
#define LILA_MATRIX_H_

#include <vector>

#include "common.h"
#include "range.h"

namespace lila {
  
  template <class coeff_t>
  class Matrix {
  public:
    using size_type = lila::size_type; 
    using coeff_type = coeff_t;    
    using vector_type = std::vector<coeff_t>;
    using iterator_t = typename vector_type::iterator;
    using const_iterator_t = typename vector_type::const_iterator;

    Matrix() : m_(0), n_(0), size_(0), data_() { };
    ~Matrix() = default;
    Matrix(const Matrix&) = default;
    Matrix& operator=(Matrix&) = default; 
    Matrix(Matrix&&) = default;
    Matrix& operator=(Matrix&&) = default;

    Matrix(size_type m, size_type n) 
      : m_(m), n_(n), size_(m*n), data_(size_) { }
    
    coeff_t operator()(size_type i, size_type j) const { return data_[i + j*m_]; }
    coeff_t& operator()(size_type i, size_type j) { return (data_[i + j*m_]); }

    void resize(size_type m, size_type n) {
      m_ = m;
      n_ = n;
      size_= m*n;
      data_.resize(size_);
    }

    void clear() {
      m_ = 0;
      n_ = 0;
      size_= 0;
      data_.clear();
    }

    range<size_type> rows() const { return range<size_type>(m_); }
    range<size_type> cols() const { return range<size_type>(n_); }

    size_type nrows() const { return m_; }
    size_type ncols() const { return n_; }
    size_type size() const { return size_; }
    
    coeff_t* data() { return data_.data(); }
    const coeff_t* data() const { return data_.data(); }

    iterator_t begin() { return data_.begin(); }
    iterator_t end() { return data_.end(); }
    const_iterator_t begin() const { return data_.begin(); }
    const_iterator_t end() const { return data_.end(); }
    const_iterator_t cbegin() const { return data_.cbegin(); }
    const_iterator_t cend() const { return data_.cend(); }
    
  private:
    size_type m_;
    size_type n_;
    size_type size_;
    vector_type data_;
  };

}

#endif
