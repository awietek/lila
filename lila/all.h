#pragma once

#if !defined(LILA_USE_MKL) && !defined(LILA_USE_ACCELERATE) && !defined(LILA_USE_LAPACK)
#error No BLAS/LAPACK backend defined for lila!
#endif

#include "common.h"
#include "vector.h"
#include "matrix.h"

#include "numeric/compare.h"
#include "numeric/complex.h"
#include "numeric/precision.h"

#include "arithmetic/matrixfunction.h"
#include "arithmetic/add.h"
#include "arithmetic/mult.h"

#include "special/random.h"
#include "special/special.h"

#include "decomp/solve.h"
#include "decomp/cholesky.h"

#include "eigen/eigen.h"
#include "eigen/eigen_sym.h"
#include "eigen/eigen_sym_tridiag.h"
#include "eigen/eigen_gen_sym.h"

#include "utils/range.h"
#include "utils/logger.h"
#include "utils/timing.h"
#include "utils/print.h"
