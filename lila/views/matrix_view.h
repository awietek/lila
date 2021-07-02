#pragma once

#include <memory>
#include <vector>

#include <lila/arithmetic/map.h>
#include <lila/arithmetic/copy.h>
#include <lila/matrix.h>
#include <lila/views/slice.h>

namespace lila {

template <class coeff_t> class Matrix;

template <class coeff_t> class MatrixView {
public:
  using vector_type = std::vector<coeff_t>;
  MatrixView() : storage_(std::make_shared<vector_type>()){};
  ~MatrixView() = default;
  MatrixView(MatrixView const&) = default;
  MatrixView(MatrixView &&) = default;
  MatrixView& operator=(MatrixView const& other) {
    printf("mv copy assing\n");
    Copy(other, *this);
    return *this;
  }
  MatrixView & operator=( MatrixView && other) {
    printf("mv move assing\n");
    Copy(other, *this);
    return *this;
  }



  MatrixView &operator=(coeff_t c) {
    Map(*this, [&c](coeff_t &x) { x = c; });
    return *this;
  }

  MatrixView(Matrix<coeff_t> const &A)
      : storage_(A.storage_), ld_(A.m()), begin_(0), m_(A.m()), n_(A.n()),
        incm_(1), incn_(1) {}

  MatrixView(Matrix<coeff_t> &A, Slice const &slice_row, Slice const &slice_col)
      : storage_(A.storage_), ld_(A.m()) {
    auto end_row = (slice_row.end == END) ? A.nrows() : slice_row.end;
    auto end_col = (slice_col.end == END) ? A.ncols() : slice_col.end;

    size_type begin = slice_row.begin + ld_ * slice_col.begin;
    size_type end = end_row + ld_ * end_col;

    assert(begin < A.size());
    assert(end < A.size());
    assert(begin <= end);

    assert(end_row >= slice_row.begin);
    assert(end_col >= slice_col.begin);

    begin_ = begin;
    m_ = (end_row - slice_row.begin) / slice_row.step;
    n_ = (end_col - slice_col.begin) / slice_col.step;
    incm_ = slice_row.step;
    incn_ = slice_col.step;
  }

  size_type m() const { return m_; }
  size_type n() const { return n_; }
  size_type size() const { return m_ * n_; }
  size_type incm() const { return incm_; }
  size_type incn() const { return incn_; }
  size_type ld() const { return ld_; }
  long use_count() const { return storage_.use_count(); }


  std::shared_ptr<vector_type> storage() { return storage_; }
  coeff_t *data() { return storage_->data() + begin_; }
  const coeff_t *data() const { return storage_->data() + begin_; }

private:
  std::shared_ptr<std::vector<coeff_t>> storage_;
  size_type ld_;
  size_type begin_;
  size_type m_, n_;
  size_type incm_, incn_;
};

} // namespace lila
