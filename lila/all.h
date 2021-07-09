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

#include "arithmetic/add.h"
#include "arithmetic/copy.h"
#include "arithmetic/dot.h"
#include "arithmetic/map.h"
#include "arithmetic/norm.h"
#include "arithmetic/scale.h"

#include "algebra/mult.h"
#include "algebra/expm.h"
#include "algebra/matrixfunction.h"

#include "special/random.h"
#include "special/special.h"
#include "special/matrix_ops.h"

#include "decomp/solve.h"
#include "decomp/cholesky.h"
#include "decomp/determinant.h"

#include "eigen/eigen.h"
#include "eigen/eigen_sym.h"
#include "eigen/eigen_sym_tridiag.h"
#include "eigen/eigen_gen_sym.h"

#include "utils/range.h"
#include "utils/logger.h"
#include "utils/timing.h"
#include "utils/print.h"
