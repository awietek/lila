#pragma once

namespace lila {

template <class coeff_t> class VectorView {
public:
  VectorView(Vector const &vec)
  {
    data_ = vec.data();
    N_ = vec.size();
    inc_ = 1;
  }

  VectorView(Vector const &vec, Slice const &slice)
  {
    assert(slice.begin < vec.size());
    assert(slice.end < vec.size());
    assert(slice.begin <= end);

    data_ = vec.data() + slice.begin;
    N_ = slice.end - slice.begin;
    inc_ = 1;
  }

  VectorView(Vector const &vec, SliceStep const &slice)
  {
    assert(slice.begin < vec.size());
    assert(slice.end < vec.size());
    assert(slice.begin <= slice.end);

    data_ = vec.data() + slice.begin;
    N_ = (slice.end - slice.begin) / slice.step;
    inc_ = step;
  }

  VectorView() = delete;
  ~VectorView() = default;
  VectorView(VectorView const &) = delete;
  VectorView &operator=(VectorView const &) = delete;
  VectorView(VectorView &&) = default;
  VectorView &operator=(VectorView &&) = default;

  inline coeff_t* data() { return data_; }
  inline int N() const { return N; }
  inline int inc() const { return inc_; }

private:
  coeff_t *data_;
  int N_;
  int inc_;
}

} // namespace lila
