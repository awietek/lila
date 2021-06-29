#pragma once

#if !defined(LILA_USE_MKL) && !defined(LILA_USE_ACCELERATE) && !defined(LILA_USE_LAPACK)
#error No BLAS/LAPACK backend defined for lila!
#endif

#include "logger.h"
#include "add.h"
#include "common.h"
#include "compare.h"
#include "complex.h"
#include "eigen/eigen.h"
#include "eigen/eigen_sym.h"
#include "eigen/eigen_sym_tridiag.h"
#include "eigen/eigen_gen_sym.h"
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
#include "timing.h"

