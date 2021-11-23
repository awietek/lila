#pragma once
#include <complex>
#include <random>

namespace lila::detail {

template <class random_engine_t, class distribution_t>
inline void random_number(random_engine_t &engine, distribution_t &distribution,
                          float &number) {
  number = distribution(engine);
}

template <class random_engine_t, class distribution_t>
inline void random_number(random_engine_t &engine, distribution_t &distribution,
                          double &number) {
  number = distribution(engine);
}

template <class random_engine_t, class distribution_t>
inline void random_number(random_engine_t &engine, distribution_t &distribution,
                          std::complex<float> &number) {
  number = {distribution(engine), distribution(engine)};
}

template <class random_engine_t, class distribution_t>
inline void random_number(random_engine_t &engine, distribution_t &distribution,
                          std::complex<double> &number) {
  number = {distribution(engine), distribution(engine)};
}

} // namespace lila::detail
