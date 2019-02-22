#include "catch.hpp"

#include <lila/allmpi.h>

template <class coeff_t>
void test_addmpi()
{
  using namespace lila;
  
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int m = 5+rank;
  int seed = 42 + 123456789*rank;

  uniform_dist_t<real_t<coeff_t>> dist(-1., 1.);
  uniform_gen_t<coeff_t> gen(dist, seed);
  
  VectorMPI<coeff_t> a(m);
  Random(a.vector_local(), gen);
  VectorMPI<coeff_t> b(m);
  Random(b.vector_local(), gen);
  
  auto res =  a + b;
  // for (int i = 0; i < size; ++i)
  //   {
  //     if (i == rank)
  // 	{
  // 	  printf("[%d] %d %d\n", rank, (int)a.size(), (int)a.size_local());
  // 	  LilaPrint(a.vector_local());
  // 	  LilaPrint(b.vector_local());
  // 	  LilaPrint(res.vector_local());
  // 	}
  //     MPI_Barrier(MPI_COMM_WORLD);
  //   }
  for (int i : a.vector_local().rows())
    REQUIRE(close(res(i), a(i) + b(i)));
  
  coeff_t fac = (coeff_t)3.1234;
  res = fac * a;
  for (int i : a.vector_local().rows())
    REQUIRE(close(res(i), fac*a(i)));

  res = a / fac;
  for (int i : a.vector_local().rows())
    REQUIRE(close(res(i), a(i) / fac));

  res = a - b;
  for (int i : a.vector_local().rows())
    REQUIRE(close(res(i), a(i) - b(i)));
  
  res = a;
  res += b; 
  for (int i : a.vector_local().rows())
    REQUIRE(close(res(i), a(i) + b(i)));

  res = a;
  res -= b; 
  for (int i : a.vector_local().rows())
    REQUIRE(close(res(i), a(i) - b(i)));

  res = a;
  res *= fac; 
  for (int i : a.vector_local().rows())
    REQUIRE(close(res(i), fac*a(i)));

  res = a;
  res /= fac; 
  for (int i : a.vector_local().rows())
    REQUIRE(close(res(i), a(i)/fac));

}

TEST_CASE( "Add MPI test", "[AddMPI]" ) {
  test_addmpi<float>();
  test_addmpi<double>();
  test_addmpi<std::complex<float>>();
  test_addmpi<std::complex<double>>();
}
