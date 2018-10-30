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

#ifndef LILA_VECTOR_H_
#define LILA_VECTOR_H_

#include <vector>

#include "common.h"
#include "range.h"

namespace lila {
  
  template <class coeff_t>
  class Vector {
  public:
    using size_type = lila::size_type; 
    using coeff_type = coeff_t;    
    using vector_type = std::vector<coeff_t>;
    using iterator_t = typename vector_type::iterator;
    using const_iterator_t = typename vector_type::const_iterator;

    Vector() : size_(0), data_() { };
    ~Vector() = default;
    Vector(const Vector&) = default;
    Vector& operator=(Vector&) = default; 
    Vector(Vector&&) = default;
    Vector& operator=(Vector&&) = default;

    explicit Vector(size_type size) 
      : size_(size), data_(size) { }
    
    coeff_t operator()(size_type i) const { return data_[i]; }
    coeff_t& operator()(size_type i) { return data_[i]; }

    void resize(size_type size) {
      size_= size;
      data_.resize(size);
    }

    void clear() {
      size_= 0;
      data_.clear();
    }

    range<size_type> rows() const { return range<size_type>(size_); }
    size_type nrows() const { return size_; }
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
    size_type size_;
    vector_type data_;
  };

}

#endif
