#pragma once

#include <memory>
#include <vector>

#include <lila/arithmetic/copy.h>
#include <lila/common.h>
#include <lila/views/slice.h>
#include <lila/views/vector_view.h>

namespace lila {

template <class coeff_t> class VectorView;

template <class coeff_t> class Vector {
public:
  using size_type = lila::size_type;
  using coeff_type = coeff_t;
  using value_type = coeff_t;
  using vector_type = std::vector<coeff_t>;
  using iterator_t = typename vector_type::iterator;
  using const_iterator_t = typename vector_type::const_iterator;

  friend class VectorView<coeff_t>;

  Vector() : storage_(std::make_shared<vector_type>()){};
  ~Vector() = default;

  Vector(Vector const &v)
      : storage_(std::make_shared<vector_type>(v.vector())){};
  Vector &operator=(Vector v) {
    swap(*this, v);
    return *this;
  }
  friend void swap(Vector &v, Vector &w) { std::swap(v.storage_, w.storage_); }
  Vector(Vector &&) = default;

  explicit Vector(size_type size)
      : storage_(std::make_shared<vector_type>(size, 0)) {}

  explicit Vector(vector_type const &vec)
      : storage_(std::make_shared<vector_type>(vec)) {}

  Vector(std::initializer_list<coeff_t> list)
      : storage_(std::make_shared<vector_type>(list.size(), 0)) {
    std::copy(list.begin(), list.end(), storage_->begin());
  }

  Vector &operator=(vector_type const &vec) {
    *storage_ = vec;
    return *this;
  };

  Vector(VectorView<coeff_t> const &view)
      : storage_(std::make_shared<vector_type>(view.size(), 0)) {
    Copy(view, VectorView<coeff_t>(*this));
  };

  Vector &operator=(VectorView<coeff_t> const &view) {
    storage_->resize(view.n());
    Copy(view, VectorView<coeff_t>(*this));
  };

  Vector &operator=(coeff_t c) {
    std::fill(storage_->begin(), storage_->end(), c);
  };

  bool operator==(Vector const &other) const {
    return *other.storage_ == *storage_;
  }

  coeff_t operator()(size_type i) const { return (*storage_)[i]; }
  coeff_t &operator()(size_type i) { return (*storage_)[i]; }

  VectorView<coeff_t> operator()(Slice const &slice) {
    return VectorView<coeff_t>(*this, slice);
  }

  operator vector_type &() { return *storage_; }

  size_type size() const { return storage_->size(); }
  size_type n() const { return storage_->size(); }
  void resize(size_type size) { storage_->resize(size, 0); }
  void clear() { storage_->clear(); }
  void push_back(coeff_t c) { storage_->push_back(c); }
  void shrink_to_fit() { storage_->shrink_to_fit(); }
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
  std::shared_ptr<vector_type> storage_;
};

} // namespace lila
