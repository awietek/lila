#pragma once

#include <lila/vector.h>
#include <lila/views/slice.h>

namespace lila {

template <class coeff_t> class Vector;

template <class coeff_t> class VectorView {
public:
  VectorView(Vector<coeff_t> &vec) {
    data_ = vec.data();
    N_ = vec.size();
    inc_ = 1;
  }

  VectorView(Vector<coeff_t> &vec, Slice && slice) {
    assert(slice.begin < vec.size());
    assert(slice.end <= vec.size());
    assert(slice.begin <= slice.end);

    data_ = vec.data() + slice.begin;
    N_ = (slice.end - slice.begin) / slice.step;
    inc_ = slice.step;
  }

  VectorView() = delete;
  ~VectorView() = default;
  VectorView(VectorView const &) = delete;
  VectorView &operator=(VectorView const &) = delete;
  VectorView(VectorView &&) = default;
  VectorView &operator=(VectorView &&) = default;

  inline coeff_t *data() { return data_; }
  inline int N() const { return N_; }
  inline int inc() const { return inc_; }

private:
  coeff_t *data_;
  int N_;
  int inc_;
};

} // namespace lila
