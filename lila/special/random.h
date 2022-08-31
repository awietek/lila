#pragma once

#include <random>

#include <lila/matrix.h>
#include <lila/numeric/complex.h>
#include <lila/vector.h>

#include <lila/decomp/solve.h>
#include <lila/detail/random_detail.h>

namespace lila {

template <class coeff_t>
using uniform_dist_t = std::uniform_real_distribution<real_t<coeff_t>>;

template <class coeff_t>
using normal_dist_t = std::normal_distribution<real_t<coeff_t>>;

template <class random_engine_t, class distribution_t, class coeff_t>
class random_generator {
public:
  using coeff_type = coeff_t;

  random_generator(distribution_t &distribution, int seed)
      : distribution_(distribution) {
    engine_.seed(seed);
  }

  coeff_t operator()() {
    coeff_t number;
    detail::random_number(engine_, distribution_, number);
    return number;
  }

  void seed(int seed) { engine_.seed(seed); }

private:
  random_engine_t engine_;
  distribution_t distribution_;
};

template <class coeff_t>
using uniform_gen_t =
    random_generator<std::mt19937, uniform_dist_t<coeff_t>, coeff_t>;
template <class coeff_t>
using normal_gen_t =
    random_generator<std::mt19937, normal_dist_t<coeff_t>, coeff_t>;

template <class coeff_t, class gen_t>
void Random(Vector<coeff_t> &vec, gen_t &gen, bool alter_generator = true) {
  if (alter_generator)
    std::for_each(vec.begin(), vec.end(), [&gen](coeff_t &c) { c = gen(); });
  else
    std::generate(vec.data(), vec.data() + vec.size(), gen);
}

template <class coeff_t, class gen_t>
void Random(Matrix<coeff_t> &mat, gen_t &gen, bool alter_generator = true) {
  if (alter_generator)
    std::for_each(mat.begin(), mat.end(), [&gen](coeff_t &c) { c = gen(); });
  else
    std::generate(mat.data(), mat.data() + mat.size(), gen);
}

template <class coeff_t> void Random(Vector<coeff_t> &vec) {
  std::random_device rd{};
  lila::normal_dist_t<coeff_t> dist(0., 1.);
  lila::normal_gen_t<coeff_t> gen(dist, rd());
  Random(vec, gen);
}

template <class coeff_t> void Random(Matrix<coeff_t> &mat) {
  std::random_device rd{};
  lila::normal_dist_t<coeff_t> dist(0., 1.);
  lila::normal_gen_t<coeff_t> gen(dist, rd());
  Random(mat, gen);
}

template <class gen_t>
Vector<typename gen_t::coeff_type> Random(lila_size_t m, gen_t &gen,
                                          bool alter_generator = true) {
  Vector<typename gen_t::coeff_type> vec(m);
  Random(vec, gen, alter_generator);
  return vec;
}

template <class gen_t>
Matrix<typename gen_t::coeff_type> Random(lila_size_t m, lila_size_t n, gen_t &gen,
                                          bool alter_generator = true) {
  Matrix<typename gen_t::coeff_type> mat(m, n);
  Random(mat, gen, alter_generator);
  return mat;
}

template <class coeff_t> Vector<coeff_t> Random(lila_size_t m) {
  Vector<coeff_t> vec(m);
  Random(vec);
  return vec;
}

template <class coeff_t> Matrix<coeff_t> Random(lila_size_t m, lila_size_t n) {
  Matrix<coeff_t> mat(m, n);
  Random(mat);
  return mat;
}

template <class matrix_t, class gen_t>
void RandomUnitary(matrix_t &mat, gen_t &gen, bool alter_generator = true) {
  assert(mat.nrows() == mat.ncols()); // Could be dropped eventually
  using coeff_t = typename matrix_t::coeff_type;
  Random(mat, gen, alter_generator);
  std::vector<coeff_t> tau = QRDecompose(mat);
  mat = QRGetQ(mat, tau);
  // TODO: get uniform Haar measure here !!!!
  // i.e. correct diagonal
  // cf. e.g. arxiv:math-ph/0609050v2
}

} // namespace lila
