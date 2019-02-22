#include "catch.hpp"
#include <mpi.h>

TEST_CASE( "MPI test", "[MPI]" ) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // printf("[%d] %d", rank, size);
  CHECK(size > 0); CHECK(rank >= 0);
}
