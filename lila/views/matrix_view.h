#pragma once

#include <lila/matrix.h>
#include <lila/views/slice.h>

namespace lila {

template <class coeff_t> class Matrix;

template <class coeff_t> class MatrixView {
public:
  MatrixView(Matrix<coeff_t> &mat)
      : ld_(mat.nrows()), M_(mat.nrows()), N_(mat.ncols()), data_(mat.data()) {}

  MatrixView(Matrix<coeff_t> &mat, Slice const &slice_row,
             Slice const &slice_col)
      : ld_(mat.nrows()) {
    auto end_row = (slice_row.end == END) ? mat.nrows() : slice_row.end; 
    auto end_col = (slice_col.end == END) ? mat.ncols() : slice_col.end; 

    size_type begin = slice_row.begin + ld_ * slice_col.begin;
    size_type end = end_row + ld_ * end_col;

    assert(begin < mat.size());
    assert(end < mat.size());
    assert(begin <= end);

    data_ = mat.data() + begin;

    assert(end_row >= slice_row.begin);
    assert(end_col >= slice_col.begin);

    M_ = end_row - slice_row.begin;
    N_ = end_col - slice_col.begin;
  }

  MatrixView() = delete;
  ~MatrixView() = default;
  MatrixView(MatrixView const &) = delete;
  MatrixView(MatrixView &&) = delete;
  MatrixView &&operator=(MatrixView &&v) = delete;
  MatrixView &&operator=(MatrixView const &v) = delete;

  inline coeff_t *data() const { return data_; }
  inline size_type M() const { return M_; }
  inline size_type N() const { return N_; }
  inline size_type ld() const { return ld_; }

private:
  size_type ld_;
  size_type M_, N_;
  coeff_t *data_;
};

} // namespace lila
