// Copyright 2018 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LILA_MPI_MPIUTILS_H_
#define LILA_MPI_MPIUTILS_H_

#include <mpi.h>

#include "../common.h"

namespace lila {

  template <class TCoeffs> inline MPI_Datatype datatype();
  template <> inline MPI_Datatype datatype<char>()
  {return MPI_CHAR;}
  template <> inline MPI_Datatype datatype<short>()
  {return MPI_SHORT;}
  template <> inline MPI_Datatype datatype<int>()
  {return MPI_INT;}
  template <> inline MPI_Datatype datatype<long>()
  {return MPI_LONG;}
  template <> inline MPI_Datatype datatype<long long>()
  {return MPI_LONG_LONG;}


  template <> inline MPI_Datatype datatype<unsigned char>()
  {return MPI_UNSIGNED_CHAR;}
  template <> inline MPI_Datatype datatype<unsigned short>()
  {return MPI_UNSIGNED_SHORT;}
  template <> inline MPI_Datatype datatype<unsigned int>()
  {return MPI_UNSIGNED;}
  template <> inline MPI_Datatype datatype<unsigned long>()
  {return MPI_UNSIGNED_LONG;}
  template <> inline MPI_Datatype datatype<unsigned long long>()
  {return MPI_UNSIGNED_LONG_LONG;}

  template <> inline MPI_Datatype datatype<float>()
  {return MPI_FLOAT;}
  template <> inline MPI_Datatype datatype<double>()
  {return MPI_DOUBLE;}

  // Allreduce
  template <class TCoeffs>
  inline int MPI_Allreduce
  (const TCoeffs *sendbuf, TCoeffs *recvbuf, int count,
   MPI_Op op, MPI_Comm comm)
  {
    MPI_Datatype type = datatype<TCoeffs>();
    return MPI_Allreduce(sendbuf, recvbuf, count, type, op, comm);
  }

  // Special implementation for complex numbers
  template <>
  inline int MPI_Allreduce<scomplex>
  (const scomplex *sendbuf, scomplex *recvbuf, int count, 
   MPI_Op op, MPI_Comm comm)
  {
    return MPI_Allreduce(sendbuf, recvbuf, count << 1, MPI_FLOAT, op, comm);
  }

  template <>
  inline int MPI_Allreduce<complex>
  (const complex *sendbuf, complex *recvbuf, int count, 
   MPI_Op op, MPI_Comm comm)
  {
    return MPI_Allreduce(sendbuf, recvbuf, count << 1, MPI_DOUBLE, op, comm);
  }

  // Alltoall
  template <class TCoeffs>
  inline int MPI_Alltoall
  (const TCoeffs *sendbuf, int sendcount, 
   TCoeffs *recvbuf, int recvcount, MPI_Comm comm)
  {
    MPI_Datatype type = datatype<TCoeffs>();
    return MPI_Alltoall(sendbuf, sendcount, type,
			 recvbuf, recvcount, type, comm);
  }
  
  template <class TCoeffs>
  inline int MPI_Alltoallv
  (const TCoeffs *sendbuf, int *sendcounts, int *sdispls,
   TCoeffs *recvbuf, int *recvcounts, int *rdispls, MPI_Comm comm)
  {
    MPI_Datatype type = datatype<TCoeffs>();
    return MPI_Alltoallv(sendbuf, sendcounts, sdispls, type,
			 recvbuf, recvcounts, rdispls, type, comm);
  }

  // Special implementation for complex numbers
  template <>
  inline int MPI_Alltoallv<scomplex>
  (const scomplex *sendbuf, int *sendcounts, int *sdispls, 
   scomplex *recvbuf, int *recvcounts, int *rdispls, MPI_Comm comm)
  {
    int n_mpi_tasks, ret;
    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_tasks);
    for (int i = 0; i < n_mpi_tasks; ++i)
      {
	sendcounts[i] <<= 1;
	sdispls[i] <<= 1;
	recvcounts[i] <<= 1;
	rdispls[i] <<= 1;
      }
    ret= MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_FLOAT, recvbuf,
		       recvcounts, rdispls, MPI_FLOAT, comm);
    for (int i = 0; i < n_mpi_tasks; ++i)
      {
	sendcounts[i] >>= 1;
	sdispls[i] >>= 1;
	recvcounts[i] >>= 1;
	rdispls[i] >>= 1;
      }
    return ret;
  }

  template <>
  inline int MPI_Alltoallv<complex>
  (const complex *sendbuf, int *sendcounts, int *sdispls, 
   complex *recvbuf, int *recvcounts, int *rdispls, MPI_Comm comm)
  {
    int n_mpi_tasks, ret;
    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_tasks);
    for (int i = 0; i < n_mpi_tasks; ++i)
      {
	sendcounts[i] <<= 1;
	sdispls[i] <<= 1;
	recvcounts[i] <<= 1;
	rdispls[i] <<= 1;
      }
    ret= MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_DOUBLE, recvbuf,
		       recvcounts, rdispls, MPI_DOUBLE, comm);
    for (int i = 0; i < n_mpi_tasks; ++i)
      {
	sendcounts[i] >>= 1;
	sdispls[i] >>= 1;
	recvcounts[i] >>= 1;
	rdispls[i] >>= 1;
      }
    return ret;
  }

  
}

#endif
