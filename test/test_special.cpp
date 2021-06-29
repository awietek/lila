#include "catch.hpp"

#include <lila/all.h>

template <class matrix_t>
bool all_zero(const matrix_t& mat)
{
  bool all_zero=true;
  for (auto i : mat.rows())
    for (auto j : mat.cols())
      all_zero &= (mat(i,j) == static_cast<typename matrix_t::coeff_type>(0.));
  return all_zero;
}

template <class matrix_t>
bool all_one(const matrix_t& mat)
{
  bool all_zero=true;
  for (auto i : mat.rows())
    for (auto j : mat.cols())
      all_zero &= (mat(i,j) == static_cast<typename matrix_t::coeff_type>(1.));
  return all_zero;
}

TEST_CASE( "special", "[core]" ) {

  int m=3;
  int n=4;
  
  // Zeros
  {
    lila::Matrix<float> fmat(m, n);
    Zeros(fmat);
    REQUIRE(all_zero(fmat));

    lila::Matrix<double> dmat(m, n);
    Zeros(dmat);
    REQUIRE(all_zero(dmat));

    lila::Matrix<std::complex<float>> cmat(m, n);
    Zeros(cmat);
    REQUIRE(all_zero(cmat));

    lila::Matrix<std::complex<double>> zmat(m, n);
    Zeros(zmat);
    REQUIRE(all_zero(zmat));
  }

  // Ones
  {
    lila::Matrix<float> fmat(m, n);
    Ones(fmat);
    REQUIRE(all_one(fmat));

    lila::Matrix<double> dmat(m, n);
    Ones(dmat);
    REQUIRE(all_one(dmat));

    lila::Matrix<std::complex<float>> cmat(m, n);
    Ones(cmat);
    REQUIRE(all_one(cmat));

    lila::Matrix<std::complex<double>> zmat(m, n);
    Ones(zmat);
    REQUIRE(all_one(zmat));
  }
  
  // Unitary
  int N = 3;
  lila::Vector<double> params(N*N);
  for (int i : lila::range<int>(N*N))
    params(i) = 0.1 * (i+1);
  auto unitary = lila::Unitary<std::complex<double>>(N, params);
  // LilaPrint(unitary);
  REQUIRE(IsUnitary(unitary));

}
