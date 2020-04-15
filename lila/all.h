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

#ifndef LILA_ALL_H_
#define LILA_ALL_H_

#if !defined(LILA_USE_MKL) && !defined(LILA_USE_ACCELERATE) && !defined(LILA_USE_LAPACK)
#error No BLAS/LAPACK backend defined for lila!
#endif

#include "logger.h"
#include "add.h"
#include "common.h"
#include "compare.h"
#include "complex.h"
#include "eigen.h"
#include "matrix.h"
#include "matrixfunction.h"
#include "mult.h"
#include "precision.h"
#include "print.h"
#include "random.h"
#include "range.h"
#include "solve.h"
#include "cholesky.h"
#include "special.h"
#include "vector.h"
#include "tmatrix.h"
#include "algorithms/expsymv.h"
#include "algorithms/lanczos.h"
#include "algorithms/bandlanczos.h"
#include "algorithms/lobpcg.h"
#include "algorithms/gramschmidt.h"
#include "timing.h"

#endif
