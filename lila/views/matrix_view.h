#pragma once

#include <memory>
#include <vector>

#include <lila/arithmetic/copy.h>
#include <lila/arithmetic/map.h>
#include <lila/matrix.h>
#include <lila/views/slice.h>

namespace lila {

template <class coeff_t> class Matrix;

template <class coeff_t> class MatrixView {
public:
  using vector_type = std::vector<coeff_t>;
  MatrixView() : storage_(std::make_shared<vector_type>()){};
  ~MatrixView() = default;
  MatrixView(MatrixView const &) = default;
  MatrixView(MatrixView &&) = default;
  MatrixView &operator=(MatrixView const &other) {
    Copy(other, *this);
    return *this;
  }
  MatrixView &operator=(MatrixView &&other) {
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
    assert(slice_row.step != 0);
    assert(slice_col.step != 0);

    auto [aslice_row, am] = adjusted_slice_length(slice_row, A.m());
    auto [aslice_col, an] = adjusted_slice_length(slice_col, A.n());

    begin_ = aslice_row.begin + ld_ * aslice_col.begin;
    m_ = am;
    n_ = an;
    incm_ = (slice_row.step > 0) ? slice_row.step : -slice_row.step;
    incn_ = (slice_col.step > 0) ? slice_col.step : -slice_col.step;
  }

  lila_size_t m() const { return m_; }
  lila_size_t n() const { return n_; }
  lila_size_t size() const { return m_ * n_; }
  lila_size_t incm() const { return incm_; }
  lila_size_t incn() const { return incn_; }
  lila_size_t ld() const { return ld_; }
  long use_count() const { return storage_.use_count(); }

  std::shared_ptr<vector_type> storage() { return storage_; }
  coeff_t *data() { return storage_->data() + begin_; }
  const coeff_t *data() const { return storage_->data() + begin_; }

private:
  std::shared_ptr<std::vector<coeff_t>> storage_;
  lila_size_t ld_;
  lila_size_t begin_;
  lila_size_t m_, n_;
  lila_size_t incm_, incn_;
};

} // namespace lila
