#pragma once

#include <memory>
#include <vector>

#include <lila/arithmetic/copy.h>
#include <lila/common.h>
#include <lila/views/matrix_view.h>
#include <lila/views/slice.h>
#include <lila/views/vector_view.h>

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

  friend class VectorView<coeff_t>;
  friend class MatrixView<coeff_t>;

  Matrix() : m_(0), n_(0), storage_(std::make_shared<vector_type>()){};
  ~Matrix() = default;

  Matrix(Matrix const &A)
      : m_(A.m_), n_(A.n_),
        storage_(std::make_shared<vector_type>(A.vector())){};

  explicit Matrix(Vector<coeff_t> const &v)
      : m_(v.size()), n_(1),
        storage_(std::make_shared<vector_type>(v.vector())){};

  Matrix &operator=(Matrix A) {
    swap(*this, A);
    return *this;
  }
  friend void swap(Matrix &A, Matrix &B) {
    using std::swap;
    swap(A.m_, B.m_);
    swap(A.n_, B.n_);
    swap(A.storage_, B.storage_);
  }
  Matrix(Matrix &&) = default;

  bool operator==(Matrix const &other) const {
    return ((other.m_ == m_) && (other.n_ == n_) &&
            (*other.storage_ == *storage_));
  }

  Matrix(size_type m, size_type n)
      : m_(m), n_(n), storage_(std::make_shared<vector_type>(m * n, 0)) {}

  Matrix(std::initializer_list<std::initializer_list<coeff_t>> listlist)
      : m_((size_type)(listlist.begin())->size()),
        n_((size_type)listlist.size()),
        storage_(std::make_shared<vector_type>(m_ * n_, 0)) {

    for (int j = 0; j < n_; j++) {
      for (int i = 0; i < m_; i++) {
        (*this)(i, j) = ((listlist.begin() + i)->begin())[j];
      }
    }
  }

  Matrix(MatrixView<coeff_t> const &view)
      : m_(view.m()), n_(view.n()),
        storage_(std::make_shared<vector_type>(m_ * n_, 0)) {
    Copy(view, MatrixView<coeff_t>(*this));
  };

  Matrix &operator=(MatrixView<coeff_t> const &A) {
    if ((m_ != A.m()) || (n_ != A.n())) {
      storage_->resize(A.m() * A.n());
      m_ = A.m();
      n_ = A.n();
    }
    Copy(A, MatrixView<coeff_t>(*this));
    return *this;
  };

  Matrix &operator=(coeff_t c) {
    std::fill(storage_->begin(), storage_->end(), c);
  };

  coeff_t operator()(size_type i, size_type j) const {
    return (*storage_)[i + j * m_];
  }
  coeff_t &operator()(size_type i, size_type j) {
    return (*storage_)[i + j * m_];
  }

  MatrixView<coeff_t> operator()(Slice const &slice_row,
                                 Slice const &slice_col) {
    return MatrixView<coeff_t>(*this, slice_row, slice_col);
  }

  VectorView<coeff_t> operator()(size_type i, Slice const &slice_col) {
    size_type end = (slice_col.end == END) ? n_ : slice_col.end;
    size_type begin = i + slice_col.begin * m_;
    size_type n = (end - slice_col.begin) / slice_col.step;
    size_type inc = slice_col.step * m_;
    return VectorView<coeff_t>(storage_, begin, n, inc);
  }

  VectorView<coeff_t> operator()(Slice const &slice_row, size_type j) {
    size_type end = (slice_row.end == END) ? m_ : slice_row.end;
    size_type begin = slice_row.begin + j * m_;
    size_type n = (end - slice_row.begin) / slice_row.step;
    size_type inc = slice_row.step;
    return VectorView<coeff_t>(storage_, begin, n, inc);
  }

  Vector<coeff_t> row(size_type i) {
    return Vector<coeff_t>(operator()(i, {0, n_}));
  }
  Vector<coeff_t> col(size_type j) {
    return Vector<coeff_t>(operator()({0, m_}, j));
  }

  size_type size() const { return storage_->size(); }
  size_type m() const { return m_; }
  size_type n() const { return n_; }
  size_type nrows() const { return m_; }
  size_type ncols() const { return n_; }
  void resize(size_type m, size_type n) {
    Matrix<coeff_t> copy = (*this);
    storage_->resize(m * n, 0);
    std::fill(storage_->begin(), storage_->end(), 0);
    for (int j = 0; j < std::min(n, n_); ++j)
      for (int i = 0; i < std::min(m, m_); ++i)
        (*storage_)[i + j * m] = copy(i, j);
    m_ = m;
    n_ = n;
  }
  void clear() {
    m_ = 0;
    n_ = 0;
    storage_->clear();
  }
  long use_count() const { return storage_.use_count(); }

  iterator_t begin() { return storage_->begin(); }
  iterator_t end() { return storage_->end(); }
  const_iterator_t begin() const { return storage_->begin(); }
  const_iterator_t end() const { return storage_->end(); }
  const_iterator_t cbegin() const { return storage_->cbegin(); }
  const_iterator_t cend() const { return storage_->cend(); }

  // Internal methods, not to be used by user
  coeff_t *data() { return storage_->data(); }
  const coeff_t *data() const { return storage_->data(); }
  std::shared_ptr<vector_type> storage() { return storage_; }
  vector_type &vector() { return *storage_; }
  vector_type const &vector() const { return *storage_; }

private:
  size_type m_;
  size_type n_;
  std::shared_ptr<vector_type> storage_;
};

} // namespace lila
