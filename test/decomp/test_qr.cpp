#include "../catch.hpp"

#include <lila/all.h>

template <class coeff_t>
void test_qr(int m, int n)
{
 
  // Test direct solve
  for (int seed : lila::range<int>(10)) {
  
    lila::uniform_dist_t<coeff_t> fdist(-1., 1.);
    lila::uniform_gen_t<coeff_t> fgen(fdist, seed);
  
    lila::Matrix<coeff_t> A(m, n);
    lila::Matrix<coeff_t> Q(m, n);
    lila::Matrix<coeff_t> R(m, n);

    Random(A, fgen);

    lila::Matrix<coeff_t> A_save(A);
    
    std::vector<coeff_t> tau = QRDecompose(A);
    Q = QRGetQ(A, tau);
    R = GetUpper(A);
    
    // Check if valid decomposition
    lila::Matrix<coeff_t> Prod(m, n);
    Mult(Q, R, Prod);
    REQUIRE(lila::close<coeff_t>(A_save, Prod));

    // check if R is upper trapezoidal
    int k = std::min(m, n);
    for (int row=0; row < k; ++row)
      for (int column=0; column < n; ++column)
	if (row > column) REQUIRE(R(row, column) == coeff_t(0.));
    
    lila::Matrix<coeff_t> Prod2(m, m);
    Mult(Q, Q, Prod2, coeff_t(1.), coeff_t(0.), 'C', 'N');

    lila::Matrix<coeff_t> Prod3(k, k);
    Mult(Q, Q, Prod3, coeff_t(1.), coeff_t(0.), 'N', 'C');

 
    lila::Matrix<coeff_t> Id2(Prod2);
    Identity(Id2);
    REQUIRE(lila::close<coeff_t>(Id2, Prod2));
    if (m <= n)
      {
	lila::Matrix<coeff_t> Id3(Prod3);
	Identity(Id3);
	REQUIRE(lila::close<coeff_t>(Id3, Prod3));
      }
  }
  
}


TEST_CASE( "qr", "[decomp]" ) {
  test_qr<float>(6,5);
  test_qr<double>(6,5);
  test_qr<std::complex<float>>(6,5);
  test_qr<std::complex<double>>(6,5);
  test_qr<float>(6,6);
  test_qr<double>(6,6);
  test_qr<std::complex<float>>(6,6);
  test_qr<std::complex<double>>(6,6);
  test_qr<float>(5,6);
  test_qr<double>(5,6);
  test_qr<std::complex<float>>(5,6);
  test_qr<std::complex<double>>(5,6);
}
