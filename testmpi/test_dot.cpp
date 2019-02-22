#include "catch.hpp"

#include <lila/allmpi.h>

TEST_CASE( "MPI Dot test", "[MPI Dot]" ) {
  using namespace lila;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int m = 5+rank;
  int seed = 42 + 123456789*rank;

  uniform_dist_t<float> sdist(-1., 1.);
  uniform_dist_t<double> ddist(-1., 1.);
  uniform_gen_t<float> sgen(sdist, seed);  
  uniform_gen_t<double> dgen(ddist, seed);
  uniform_gen_t<scomplex> cgen(sdist, seed);  
  uniform_gen_t<complex> zgen(ddist, seed);

  VectorMPI<complex> zvec(m);
  Random(zvec.vector_local(), zgen);
  VectorMPI<complex> zvec2(m);
  Random(zvec2.vector_local(), zgen);
  // for (int i = 0; i < size; ++i)
  //   if (i == rank)
  //     {
  // 	printf("[%d] %d %d\n", rank, (int)vec.size(), (int)vec.size_local());
  // 	LilaPrint(vec.vector_local());
  // 	LilaPrint(vec2.vector_local());
  // 	MPI_Barrier(MPI_COMM_WORLD);
  //     }
  complex zdot = Dot(zvec, zvec2);

  complex zdot_test = Dot(zvec.vector_local(), zvec2.vector_local());  
  double zdot_test_real = real(zdot_test);
  double zdot_test_imag = imag(zdot_test);
  MPI_Allreduce(MPI_IN_PLACE, &zdot_test_real, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &zdot_test_imag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  zdot_test = complex(zdot_test_real, zdot_test_imag);
  REQUIRE(close(zdot, zdot_test));

  // LilaPrint(zdot);
  // LilaPrint(zdot_test);

  VectorMPI<float> svec(m);
  Random(svec.vector_local(), sgen);
  VectorMPI<float> svec2(m);
  Random(svec2.vector_local(), sgen);
  float sdot = Dot(svec, svec2);
  float sdot_test = Dot(svec.vector_local(), svec2.vector_local());  
  MPI_Allreduce(MPI_IN_PLACE, &sdot_test, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  REQUIRE(close(sdot, sdot_test));

  // LilaPrint(sdot);
  // LilaPrint(sdot_test);

  VectorMPI<double> dvec(m);
  Random(dvec.vector_local(), dgen);
  VectorMPI<double> dvec2(m);
  Random(dvec2.vector_local(), dgen);
  double ddot = Dot(dvec, dvec2);
  double ddot_test = Dot(dvec.vector_local(), dvec2.vector_local());  
  MPI_Allreduce(MPI_IN_PLACE, &ddot_test, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  REQUIRE(close(ddot, ddot_test));

  // LilaPrint(ddot);
  // LilaPrint(ddot_test);

  VectorMPI<scomplex> cvec(m);
  Random(cvec.vector_local(), cgen);
  VectorMPI<scomplex> cvec2(m);
  Random(cvec2.vector_local(), cgen);
  scomplex cdot = Dot(cvec, cvec2);
  scomplex cdot_test = Dot(cvec.vector_local(), cvec2.vector_local());  
  float cdot_test_real = real(cdot_test);
  float cdot_test_imag = imag(cdot_test);
  MPI_Allreduce(MPI_IN_PLACE, &cdot_test_real, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &cdot_test_imag, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  cdot_test = scomplex(cdot_test_real, cdot_test_imag);
  REQUIRE(close(cdot, cdot_test));

  // LilaPrint(cdot);
  // LilaPrint(cdot_test);
  
  float snorm = Norm(svec);
  float snorm_test = pow(Norm(svec.vector_local()), 2);
  MPI_Allreduce(MPI_IN_PLACE, &snorm_test, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  snorm_test = sqrt(snorm_test);
  REQUIRE(close(snorm, snorm_test));
  // LilaPrint(snorm);
  // LilaPrint(snorm_test);

  double dnorm = Norm(dvec);
  double dnorm_test = pow(Norm(dvec.vector_local()), 2);
  MPI_Allreduce(MPI_IN_PLACE, &dnorm_test, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  dnorm_test = sqrt(dnorm_test);
  REQUIRE(close(dnorm, dnorm_test));
  // LilaPrint(dnorm);
  // LilaPrint(dnorm_test);

  float cnorm = Norm(cvec);
  float cnorm_test = pow(Norm(cvec.vector_local()), 2);
  MPI_Allreduce(MPI_IN_PLACE, &cnorm_test, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  cnorm_test = sqrt(cnorm_test);
  REQUIRE(close(cnorm, cnorm_test));
  // LilaPrint(cnorm);
  // LilaPrint(cnorm_test);

  double znorm = Norm(zvec);
  double znorm_test = pow(Norm(zvec.vector_local()), 2);
  MPI_Allreduce(MPI_IN_PLACE, &znorm_test, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  znorm_test = sqrt(znorm_test);
  REQUIRE(close(znorm, znorm_test));
  // LilaPrint(znorm);
  // LilaPrint(znorm_test);

}
