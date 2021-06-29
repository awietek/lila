#pragma once

namespace lila {

template <class coeff_t> class MatrixView {
public:
  MatrixView(Matrix const &mat)
      : ld_(mat.nrows()), M_(mat.nrows()), N_(mat.ncols()), data_(mat.data()) {}

  MatrixView(Matrix const &mat, Slice const &slice_row, Slice const &slice_col)
      : ld_(mat.nrows()) {
    int begin = slice_row.begin + ld_ * slice_col.begin;
    int end = slice_row.end + ld_ * slice_col.end;

    assert(begin < mat.size());
    assert(end < mat.size());
    assert(begin <= end);

    data_ = mat.data() + begin;

    assert(slice_row.end >= slice_row.begin);
    assert(slice_col.end >= slice_col.begin);

    M_ = slice_row.end - slice_row.begin;
    N_ = slice_col.end - slice_col.begin;
  }

  MatrixView() = delete;
  ~MatrixView() = default;
  MatrixView(MatrixView const &) = delete;
  MatrixView &operator=(MatrixView const &) = delete;
  MatrixView(MatrixView &&) = default;
  MatrixView &operator=(MatrixView &&) = default;

  inline coeff_t *data() { return data_; }
  inline int M() const { return M_; }
  inline int N() const { return N_; }
  inline int ld() const { return ld_; }

private:
  coeff_t *data_;
  int M_, N_;
  int ld_;
}

} // namespace lila
