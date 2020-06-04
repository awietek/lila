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
#include <fstream>
#include <sstream>
#include <iterator>
#include <iostream>
#include <iomanip>

#include "common.h"
#include "range.h"
#include "strings.h"

namespace lila {
  
  template <class coeff_t>
  class Vector {
  public:
    using size_type = lila::size_type; 
    using coeff_type = coeff_t; 
    using value_type = coeff_t;   
    using vector_type = std::vector<coeff_t>;
    using iterator_t = typename vector_type::iterator;
    using const_iterator_t = typename vector_type::const_iterator;

    Vector() : size_(0), data_() { };
    ~Vector() = default;
    Vector(const Vector&) = default;
    Vector& operator=(const Vector&) = default; 
    Vector(Vector&&) = default;
    Vector& operator=(Vector&&) = default;


    explicit Vector(size_type size) 
      : size_(size), data_(size, 0) { }
    explicit Vector(const vector_type& vec) 
      : size_(vec.size()), data_(vec) { }
    Vector& operator=(std::vector<coeff_t>& vec)
    {
      data_ = vec;
      size_ = (int)vec.size();
      return *this;
    };

    bool operator==(Vector const& other)
    { return (other.data_ == data_) && (other.size_ == size_); }

    
    coeff_t operator()(size_type i) const { return data_[i]; }
    coeff_t& operator()(size_type i) { return data_[i]; }

    operator std::vector<coeff_t>& () { return data_; }

    void resize(size_type size) {
      size_= size;
      data_.resize(size);
    }

    void clear() {
      size_= 0;
      data_.clear();
    }
    void shrink_to_fit() { data_.shrink_to_fit(); }

    range<size_type> rows() const { return range<size_type>(size_); }
    size_type nrows() const { return size_; }
    int nblocks() const { return 1; }    
    size_type size() const { return size_; }
    
    coeff_t* data() { return data_.data(); }
    const coeff_t* data() const { return data_.data(); }

    iterator_t begin() { return data_.begin(); }
    iterator_t end() { return data_.end(); }
    const_iterator_t begin() const { return data_.begin(); }
    const_iterator_t end() const { return data_.end(); }
    const_iterator_t cbegin() const { return data_.cbegin(); }
    const_iterator_t cend() const { return data_.cend(); }

    void push_back(const coeff_t& c) { 
      ++size_;
      data_.push_back(c);
    }
    
  private:
    size_type size_;
    vector_type data_;
  };

  template <class coeff_t>
  Vector<coeff_t> String2Vector(const std::string& str)
  {
    std::istringstream stream(str);
    std::vector<std::string> split(std::istream_iterator<std::string>{stream},
				   std::istream_iterator<std::string>());

    Vector<coeff_t> vector(split.size());
    int i=0;
    for (auto str : split)
      vector(i++) = string2number<coeff_t>(str);
    return vector;
  }

  template <class coeff_t>
  std::string Vector2String(const Vector<coeff_t>& vector, bool vertical = true)
  {
    char whitespace = vertical ? '\n' : ' ';
    std::stringstream ss;
    ss << std::setprecision(16);
    for (auto i : vector.rows())
      ss << vector(i) << whitespace;
    return ss.str();
  }

  template <class coeff_t>
  Vector<coeff_t> ReadVector(const std::string& filename)
  {
    std::ifstream t(filename);

    if(t.fail()) 
      {
	std::cerr << "Lila, Error in ReadVector: " 
		  << "Could not open file with filename ["
		  << filename << "] given. Abort." << std::endl;
	exit(EXIT_FAILURE);
      }

    std::string str((std::istreambuf_iterator<char>(t)),
		    std::istreambuf_iterator<char>());
    return String2Vector<coeff_t>(str);
  }

  template <class coeff_t>
  void WriteVector(const Vector<coeff_t>& vector, std::string filename, 
		   bool vertical = true)
  {
    std::ofstream t;
    t.open(filename);

    if(t.fail()) 
      {
	std::cerr << "Lila, Error in WriteVector: " 
		  << "Could not open file with filename ["
		  << filename << "] given. Abort." << std::endl;
	exit(EXIT_FAILURE);
      }

    t << Vector2String(vector, vertical);
    t.close();
  }

}

#endif
