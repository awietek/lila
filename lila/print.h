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

#ifndef LILA_PRINT_H_
#define LILA_PRINT_H_

#include <complex>

#include "matrix.h"
#include "vector.h"

#define LilaPrint(X) lila::PrintPretty(#X,X)

namespace lila {

  // TODO: allow for different precision outputs
  inline void PrintPretty(const char* identifier, const Matrix<float>& mat) {
    printf("%s:\n", identifier);
    for (auto i : mat.rows())
      {
	for (auto j : mat.cols())
	  printf("%10.8g ", mat(i,j));
	printf("\n");
      }
    printf("\n");
  }

  inline void PrintPretty(const char* identifier, const Matrix<double>& mat) {
    printf("%s:\n", identifier);
    for (auto i : mat.rows())
      {
	for (auto j : mat.cols())
	  printf("%10.8g ", mat(i,j));
	printf("\n");
      }
    printf("\n");
  }

  inline void PrintPretty(const char* identifier, 
			  const Matrix<std::complex<float>>& mat) {
    printf("%s:\n", identifier);
    for (auto i : mat.rows())
      {
	for (auto j : mat.cols())
	  printf("%10.8g%-+8.8gj ", mat(i,j).real(), mat(i,j).imag());
	printf("\n");
      }
    printf("\n");
  }


  inline void PrintPretty(const char* identifier, 
			  const Matrix<std::complex<double>>& mat) {
    printf("%s:\n", identifier);
    for (auto i : mat.rows())
      {
	for (auto j : mat.cols())
	  printf("%10.8g%-+8.8gj ", mat(i,j).real(), mat(i,j).imag());
	printf("\n");
      }
    printf("\n");
  }

  
  inline void PrintPretty(const char* identifier, const Vector<float>& vec) {
    printf("%s:\n", identifier);
    for (auto i : vec.rows())
      printf("%10.8g ", vec(i));
    printf("\n");
  }

  inline void PrintPretty(const char* identifier, const Vector<double>& vec) {
    printf("%s:\n", identifier);
    for (auto i : vec.rows())
      printf("%10.8g ", vec(i));
    printf("\n");
  }

  inline void PrintPretty(const char* identifier, 
  			  const Vector<std::complex<float>>& vec) {
    printf("%s:\n", identifier);
    for (auto i : vec.rows())
      printf("%10.8g%-+8.8gj ", vec(i).real(), vec(i).imag());
    printf("\n");
  }


  inline void PrintPretty(const char* identifier, 
  			  const Vector<std::complex<double>>& vec) {
    printf("%s:\n", identifier);
    for (auto i : vec.rows())
      printf("%10.8g%-+8.8gj ", vec(i).real(), vec(i).imag());
    printf("\n");
  }
  
}

#endif
