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

#include <cassert>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include <vector>
#include <iterator>

#include "vector.h"
#include "common.h"
#include "range.h"

namespace lila {
  
  template <class coeff_t>
  class Matrix {
  public:
    using size_type = lila::size_type; 
    using coeff_type = coeff_t;
    using value_type = coeff_t;        
    using vector_type = std::vector<coeff_t>;
    using iterator_t = typename vector_type::iterator;
    using const_iterator_t = typename vector_type::const_iterator;

    Matrix() : m_(0), n_(0), size_(0), data_() { };
    ~Matrix() = default;
    Matrix(const Matrix&) = default;
    Matrix& operator=(const Matrix&) = default; 
    Matrix(Matrix&&) = default;
    Matrix& operator=(Matrix&&) = default;

   bool operator==(Matrix const& other) const
    { return ((other.data_ == data_) && (other.size_ == size_) &&
	      (other.m_ == m_) && (other.n_ == n_)); }
    
    Matrix(size_type m, size_type n) 
      : m_(m), n_(n), size_(m*n), data_(size_, 0) { }

    coeff_t operator()(size_type i, size_type j) const
    { return data_[i + j*m_]; }
    coeff_t& operator()(size_type i, size_type j)
    { return (data_[i + j*m_]); }

    void resize(size_type m, size_type n) {
      std::vector<coeff_t> data_copy = data_;
      size_= m*n;
      data_.resize(size_);
      std::fill(data_.begin(), data_.end(), 0.);

      for (int i=0; i < std::min(m_, m); ++i)
      	for (int j=0; j < std::min(n_, n); ++j)
      	  data_[i + j*m] = data_copy[i + j*m_];

      m_ = m;
      n_ = n;
    }

    void clear() {
      m_ = 0;
      n_ = 0;
      size_= 0;
      data_.clear();
    }

    Vector<coeff_t> row(const int& i) const
    {
      Vector<coeff_t> rowi(n_);
      for (auto j : cols())
	rowi(j) = (*this)(i, j);
      return rowi;
    }
    Vector<coeff_t> col(const int& j) const
    {
      Vector<coeff_t> coli(n_);
      for (auto i : rows())
	coli(i) = (*this)(i, j);
      return coli;
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

  template <class coeff_t>
  coeff_t Trace(const Matrix<coeff_t>& mat) 
  {
    // TODO: optimize ??
    coeff_t res = 0;
    for (size_type i : range<size_type>(std::min(mat.nrows(), mat.ncols())) )
      res += mat(i,i);
    return res;
  }
  
  template <class coeff_t>
  Matrix<coeff_t> Transpose(const Matrix<coeff_t>& mat) 
  {
    Matrix<coeff_t> mat_t(mat.ncols(), mat.nrows());
    for (auto i : mat.rows())
      for (auto j : mat.cols())
	mat_t(j, i) = mat(i, j);
    return mat_t;
  }

  template <class coeff_t>
  Vector<coeff_t> Diag(const Matrix<coeff_t>& mat) 
  {
    long size = std::min(mat.nrows(), mat.ncols());
    Vector<coeff_t> vec_t(size);
    for (long i=0; i<size; ++i)
	vec_t(i) = mat(i, i);
    return vec_t;
  }



  template <class coeff_t>
  Matrix<coeff_t> ParseMatrix(const std::string& str)
  {
    using size_type = lila::size_type; 

    std::istringstream stream(str);
    std::vector<std::string> split(std::istream_iterator<std::string>{stream},
				   std::istream_iterator<std::string>());
    std::string dimstring = split[0];
    unsigned open = dimstring.find('[');
    unsigned comma = dimstring.find(',');   
    unsigned close = dimstring.find(']');  
    
    size_type m = std::stoi(dimstring.substr(open+1, comma - open));
    size_type n = std::stoi(dimstring.substr(comma+1, close - comma));    
    
    unsigned dim = static_cast<unsigned>(m*n);
    assert(split.size() == dim + 1);

    Matrix<coeff_t> matrix(m,n);
    for (auto i : range<size_type>(m))
      for (auto j : range<size_type>(n))
	{
	  std::istringstream is(split[i*m + j + 1]);
	  is >> matrix(i,j);
	}      
    return matrix;
  }

  template <class coeff_t>
  std::string WriteMatrix(const Matrix<coeff_t>& matrix)
  {
    std::stringstream ss;
    ss << " [" << matrix.nrows() << "," << matrix.ncols() << "] "; 
    ss << std::setprecision(18);
    for (auto i : matrix.rows())
      for (auto j : matrix.cols())
	ss << matrix(i,j) << " ";
    return ss.str();
  }

}

#endif
