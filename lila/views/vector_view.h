#pragma once

#include <lila/arithmetic/copy.h>
#include <lila/vector.h>
#include <lila/views/slice.h>

namespace lila {

template <class coeff_t> class Vector;
template <class coeff_t> class VectorView {
public:
  using vector_type = std::vector<coeff_t>;

  VectorView() : storage_(std::make_shared<vector_type>()){};
  ~VectorView() = default;
  VectorView(VectorView const &) = default;
  VectorView(VectorView &&) = default;
  VectorView &operator=(VectorView const &other) {
    Copy(other, *this);
    return *this;
  }
  VectorView &operator=(VectorView &&other) {
    Copy(other, *this);
    return *this;
  }

  VectorView &operator=(coeff_t c) {
    Map(*this, [&c](coeff_t &x) { x = c; });
    return *this;
  }

  VectorView(Vector<coeff_t> const &v)
      : storage_(v.storage_), begin_(0), n_(v.size()), inc_(1) {}

  VectorView(Vector<coeff_t> &v, Slice const &slice) : storage_(v.storage_) {

    assert(slice.step != 0);
    size_type length = v.size();
    auto [aslice, alength] = adjusted_slice_length(slice, length);

    begin_ = aslice.begin;
    n_ = alength;
    inc_ = (slice.step > 0) ? slice.step : -slice.step;


    assert(aslice.begin < v.size());
    assert(aslice.end <= v.size());
    assert(aslice.begin <= aslice.end);

    // printf("slice: b: %ld, e: %ld, s: %ld, l: %ld\n", slice.begin, slice.end,
    //        slice.step, slice.end - slice.begin);
    // printf("alice: b: %ld, e: %ld, s: %ld, l: %ld\n", aslice.begin, aslice.end,
    //        aslice.step, alength);
    // printf("b: %ld, n: %ld, i: %ld\n\n", begin_, n_, inc_);
  }

  VectorView(std::shared_ptr<vector_type> const &storage, size_type begin,
             size_type n, size_type inc)
      : storage_(storage), begin_(begin), n_(n), inc_(inc) {}

  size_type size() const { return n_; }
  size_type n() const { return n_; }
  size_type inc() const { return inc_; }
  long use_count() const { return storage_.use_count(); }

  std::shared_ptr<vector_type> storage() { return storage_; }
  coeff_t *data() { return storage_->data() + begin_; }
  const coeff_t *data() const { return storage_->data() + begin_; }

private:
  std::shared_ptr<vector_type> storage_;
  size_type begin_;
  size_type n_;
  size_type inc_;
};

} // namespace lila
