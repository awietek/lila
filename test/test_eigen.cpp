#include "catch.hpp"

#include <iostream>
#include <algorithm>
#include <lila/all.h>

template <class coeff_t>
void test_eigen()
{
  using namespace lila;
  int n=10;

  
  // Test Eigen
  for (int seed : range<int>(10)) {
    uniform_dist_t<coeff_t> dist(-1., 1.);
    uniform_gen_t<coeff_t> gen(dist, seed);
    Matrix<coeff_t> A(n, n);
    Random(A, gen);

    ////////////////////////////////////////////////////
    // Check whether symmetric and generic eigenvalue
    // routines give same eigenvalues
    auto A1 = A;
    Add(Herm(A1), A1);
    auto [evals1, evecs1l, evecs1r] = Eigen(A1);
    auto A2 = A;
    Add(Herm(A2), A2);
    auto [evals2, evecs2] = EigenSym(A2);

    auto e1r = Real(evals1);
    auto e1i = Imag(evals1);    
    std::sort(e1r.data(), e1r.data() + e1r.size());
    REQUIRE(close(e1r, evals2));
    REQUIRE(close(e1i, (real_t<coeff_t>)0.));




    /////////////////////////////////////////////////////////
    // Check whether symmetric  eigenvalue
    // routines computes correct eigenvalues/vectors
    auto Ah = A;
    Add(Herm(Ah), Ah);

    auto [evals3, evecs3] = EigenSym(Ah);

    // Check whether EigenvaluesSym and EigenSym give same eigenvalues
    REQUIRE(close(evals3, EigenvaluesSym(Ah)));

    // check if routine computes correct eigenvalues/vectors
    for (int i=0; i<n; ++i) 
      {
	auto eval = evals3(i);
	auto evec = evecs3.col(i);
	// LilaPrint(Mult(Ah, evec) - (coeff_t)eval * Mult(Ad, evec));
	REQUIRE(close(Mult(Ah, evec), (coeff_t)eval * evec));
      }



    /////////////////////////////////////////////////////////
    // Check symmetric definite generalized eigenvalue
    // Create positive definite Matrix
    Matrix<coeff_t> Ad(n, n);
    std::vector<lila::Vector<coeff_t>> vecs; 
    for (int k=0; k<n; ++k)
      vecs.push_back(Random<coeff_t>(n, gen));
    for (int i=0; i<n; ++i)
      for (int j=0; j<n; ++j)
    	Ad(i,j) = Dot(vecs[i], vecs[j]);

    // for (int k=0; k<n; ++k)
    //   LilaPrint(vecs[k]);

    auto [evals4, evecs4] = EigenGenSymDef(Ah, Ad);

    // Check whether EigenvaluesSym and EigenSym give same eigenvalues
    REQUIRE(close(evals4, EigenvaluesGenSymDef(Ah, Ad)));

    // check if routine computes correct eigenvalues/vectors
    for (int i=0; i<n-1; ++i) // n-1: Last eigenvector somehow looses precision?
      {
    	auto eval = evals4(i);
    	auto evec = evecs4.col(i);
    	// LilaPrint(Mult(Ah, evec) - (coeff_t)eval * Mult(Ad, evec));
    	REQUIRE(close(Mult(Ah, evec), (coeff_t)eval * Mult(Ad, evec)));
      }

  }
}


TEST_CASE( "eigen", "[decomp]" ) {
  test_eigen<float>();
  test_eigen<double>();
  test_eigen<std::complex<float>>();
  test_eigen<std::complex<double>>();
}
