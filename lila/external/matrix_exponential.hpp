#ifndef EXPOKIT_MATRIX_EXPONENTIAL_HPP_
#define EXPOKIT_MATRIX_EXPONENTIAL_HPP_

#include <complex>

using namespace std;
//
//  Complex functions.
//
inline complex <double> *c8mat_expm1 ( int n, complex <double> a[] );
//
//  Real functions.
//
inline double *r8mat_expm1 ( int n, double a[] );
inline double *r8mat_expm2 ( int n, double a[] );
inline double *r8mat_expm3 ( int n, double a[] );

#include "matrix_exponential.ih"

#endif
