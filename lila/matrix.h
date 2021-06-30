#pragma once

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <sstream>

#include <iterator>
#include <vector>

#include <lila/common.h>
#include <lila/vector.h>
#include <lila/views/matrix_view.h>
#include <lila/views/slice.h>

namespace lila {

template <class coeff_t> class Vector;
template <class coeff_t> class VectorView;
template <class coeff_t> class MatrixView;

template <class coeff_t> class Matrix {
public:
  using size_type = lila::size_type;
  using coeff_type = coeff_t;
  using value_type = coeff_t;
  using vector_type = std::vector<coeff_t>;
  using iterator_t = typename vector_type::iterator;
  using const_iterator_t = typename vector_type::const_iterator;

  Matrix() : m_(0), n_(0), size_(0), data_(){};
  ~Matrix() = default;
  Matrix(const Matrix &) = default;
  Matrix &operator=(const Matrix &) = default;
  Matrix(Matrix &&) = default;
  Matrix &operator=(Matrix &&) = default;

  bool operator==(Matrix const &other) const {
    return ((other.data_ == data_) && (other.size_ == size_) &&
            (other.m_ == m_) && (other.n_ == n_));
  }

  Matrix(size_type m, size_type n)
      : m_(m), n_(n), size_(m * n), data_(size_, 0) {}

  Matrix(std::initializer_list<std::initializer_list<coeff_t>> listlist)
      : m_((size_type)(listlist.begin())->size()),
        n_((size_type)listlist.size()), size_(m_ * n_), data_(size_, 0) {
    for (int i = 0; i < m_; i++) {
      for (int j = 0; j < n_; j++) {
        data_[i + j * m_] = ((listlist.begin() + i)->begin())[j];
      }
    }
  }

  Matrix(MatrixView<coeff_t> const &view)
      : m_(view.M()), n_(view.N()), size_(m_ * n_), data_(size_) {
    Copy(std::move(view), MatrixView<coeff_t>(*this));
  };
  Matrix &operator=(MatrixView<coeff_t> const &view) {
    size_ = view.M() * view.N();
    data_.resize(size_);
    Copy(view, *this);
  };

  coeff_t operator()(size_type i, size_type j) const {
    return data_[i + j * m_];
  }
  coeff_t &operator()(size_type i, size_type j) { return (data_[i + j * m_]); }

  MatrixView<coeff_t> operator()(Slice const &slice_row, Slice const &slice_col) {
    return MatrixView<coeff_t>(*this, slice_row, slice_col);
  }

  VectorView<coeff_t> operator()(size_type i, Slice const &slice_col) {
    auto end = (slice_col.end == END) ? n_ : slice_col.end; 
    coeff_t *vdata = data_.data() + i + slice_col.begin * m_;
    size_type vN = (end - slice_col.begin) / slice_col.step;
    size_type vinc = slice_col.step * m_;
    return VectorView<coeff_t>{vdata, vN, vinc};
  }

  VectorView<coeff_t> operator()(Slice const &slice_row, size_type j) {
    auto end = (slice_row.end == END) ? m_ : slice_row.end; 
    coeff_t *vdata = data_.data() + slice_row.begin + j * m_;
    size_type vN = (end - slice_row.begin) / slice_row.step;
    size_type vinc = slice_row.step;
    return VectorView<coeff_t>{vdata, vN, vinc};
  }


  void resize(size_type m, size_type n) {
    std::vector<coeff_t> data_copy = data_;
    size_ = m * n;
    data_.resize(size_);
    std::fill(data_.begin(), data_.end(), 0.);

    for (int i = 0; i < std::min(m_, m); ++i)
      for (int j = 0; j < std::min(n_, n); ++j)
        data_[i + j * m] = data_copy[i + j * m_];

    m_ = m;
    n_ = n;
  }

  void clear() {
    m_ = 0;
    n_ = 0;
    size_ = 0;
    data_.clear();
  }

  Vector<coeff_t> row(const int &i) const {
    Vector<coeff_t> rowi(n_);
    for (size_type j = 0; j < n_; ++j)
      rowi(j) = (*this)(i, j);
    return rowi;
  }
  Vector<coeff_t> col(const int &j) const {
    Vector<coeff_t> coli(n_);
    for (size_type i = 0; i < m_; ++i)
      coli(i) = (*this)(i, j);
    return coli;
  }

  size_type nrows() const { return m_; }
  size_type ncols() const { return n_; }
  size_type size() const { return size_; }

  coeff_t *data() { return data_.data(); }
  const coeff_t *data() const { return data_.data(); }

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

template <class coeff_t> coeff_t Trace(Matrix<coeff_t> const &mat) {
  // TODO: optimize ??
  coeff_t res = 0;
  for (size_type i = 0; i < std::min(mat.nrows(), mat.ncols()); ++i)
    res += mat(i, i);
  return res;
}

template <class coeff_t> Matrix<coeff_t> Transpose(Matrix<coeff_t> const &mat) {
  Matrix<coeff_t> mat_t(mat.ncols(), mat.nrows());
  for (int i = 0; i < mat.nrows(); ++i)
    for (int j = 0; j < mat.ncols(); ++j)
      mat_t(j, i) = mat(i, j);
  return mat_t;
}

template <class coeff_t> Vector<coeff_t> Diag(Matrix<coeff_t> const &mat) {
  long size = std::min(mat.nrows(), mat.ncols());
  Vector<coeff_t> vec_t(size);
  for (long i = 0; i < size; ++i)
    vec_t(i) = mat(i, i);
  return vec_t;
}

template <class coeff_t> Matrix<coeff_t> ParseMatrix(std::string const &str) {
  using size_type = lila::size_type;

  std::istringstream stream(str);
  std::vector<std::string> split(std::istream_iterator<std::string>{stream},
                                 std::istream_iterator<std::string>());
  std::string dimstring = split[0];
  unsigned open = dimstring.find('[');
  unsigned comma = dimstring.find(',');
  unsigned close = dimstring.find(']');

  size_type m = std::stoi(dimstring.substr(open + 1, comma - open));
  size_type n = std::stoi(dimstring.substr(comma + 1, close - comma));

  // unsigned dim = static_cast<unsigned>(m * n);
  // assert(split.size() == dim + 1);

  Matrix<coeff_t> matrix(m, n);
  for (size_type i = 0; i < m; ++i)
    for (size_type j = 0; j < n; ++j) {
      std::istringstream is(split[i * m + j + 1]);
      is >> matrix(i, j);
    }
  return matrix;
}

template <class coeff_t>
std::string WriteMatrix(Matrix<coeff_t> const &matrix) {
  std::stringstream ss;
  ss << " [" << matrix.nrows() << "," << matrix.ncols() << "] ";
  ss << std::setprecision(18);
  for (auto i : matrix.rows())
    for (auto j : matrix.cols())
      ss << matrix(i, j) << " ";
  return ss.str();
}

} // namespace lila
