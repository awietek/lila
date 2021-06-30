#pragma once

#include <lila/arithmetic/copy.h>
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

  VectorView(Vector<coeff_t> &vec, Slice const &slice) {
    auto end = (slice.end == END) ? vec.size() : slice.end; 
    assert(slice.begin < vec.size());
    assert(end <= vec.size());
    assert(slice.begin <= slice.end);

    data_ = vec.data() + slice.begin;
    N_ = (end - slice.begin) / slice.step;
    inc_ = slice.step;
  }

  VectorView(coeff_t *data, size_type N, size_type inc)
      : data_(data), N_(N), inc_(inc) {}

  VectorView() = delete;
  ~VectorView() = default;
  VectorView(VectorView const &) = delete;
  VectorView(VectorView &&) = delete;
  VectorView &operator=(VectorView const &v) = delete;
  VectorView &operator=(VectorView &&v) = delete;

  inline coeff_t *data() const { return data_; }
  inline size_type N() const { return N_; }
  inline size_type inc() const { return inc_; }

private:
  coeff_t *data_;
  size_type N_;
  size_type inc_;
};

} // namespace lila
