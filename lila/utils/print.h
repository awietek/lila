#pragma once

#include <complex>

#include <lila/common.h>
#include <lila/matrix.h>
#include <lila/vector.h>

#define LilaPrint(X) lila::PrintPretty(#X, X)

namespace lila {

// TODO: allow for different precision outputs
inline void PrintPretty(const char *identifier, const float &number) {
  printf("%s:\n", identifier);
  printf("%10.8g\n", number);
}

inline void PrintPretty(const char *identifier, const double &number) {
  printf("%s:\n", identifier);
  printf("%18.16g\n", number);
}

inline void PrintPretty(const char *identifier,
                        const std::complex<float> &number) {
  printf("%s:\n", identifier);
  printf("%10.8g%-+8.8gj\n", number.real(), number.imag());
}

inline void PrintPretty(const char *identifier,
                        const std::complex<double> &number) {
  printf("%s:\n", identifier);
  printf("%18.16g%-+18.16gj\n", number.real(), number.imag());
}

inline void PrintPretty(const char *identifier, const int32_t &number) {
  printf("%s:\n", identifier);
  printf("%d\n", number);
}

inline void PrintPretty(const char *identifier, const uint32_t &number) {
  printf("%s:\n", identifier);
  printf("%d\n", number);
}

inline void PrintPretty(const char *identifier, const int64_t &number) {
  printf("%s:\n", identifier);
  printf("%lld\n", number);
}

inline void PrintPretty(const char *identifier, const uint64_t &number) {
  printf("%s:\n", identifier);
  printf("%llu\n", number);
}

inline void PrintPretty(const char *identifier, const Matrix<float> &mat) {
  printf("%s:\n", identifier);
  for (int i = 0; i < mat.nrows(); ++i) {
    for (int j = 0; j < mat.ncols(); ++j)
      printf("%10.8g ", mat(i, j));
    printf("\n");
  }
  printf("\n");
}

inline void PrintPretty(const char *identifier, const Matrix<double> &mat) {
  printf("%s:\n", identifier);
  for (int i = 0; i < mat.nrows(); ++i) {
    for (int j = 0; j < mat.ncols(); ++j)
      printf("%10.8g ", mat(i, j));
    printf("\n");
  }
  printf("\n");
}

inline void PrintPretty(const char *identifier,
                        const Matrix<std::complex<float>> &mat) {
  printf("%s:\n", identifier);
  for (int i = 0; i < mat.nrows(); ++i) {
    for (int j = 0; j < mat.ncols(); ++j)
      printf("%10.8g%-+8.8gj ", mat(i, j).real(), mat(i, j).imag());
    printf("\n");
  }
  printf("\n");
}

inline void PrintPretty(const char *identifier,
                        const Matrix<std::complex<double>> &mat) {
  printf("%s:\n", identifier);
  for (int i = 0; i < mat.nrows(); ++i) {
    for (int j = 0; j < mat.ncols(); ++j)
      printf("%10.8g%-+8.8gj ", mat(i, j).real(), mat(i, j).imag());
    printf("\n");
  }
  printf("\n");
}

inline void PrintPretty(const char *identifier, const Vector<float> &vec) {
  printf("%s:\n", identifier);
  for (int i = 0; i < vec.size(); ++i)
    printf("%10.8g ", vec(i));
  printf("\n");
}

inline void PrintPretty(const char *identifier, const Vector<double> &vec) {
  printf("%s:\n", identifier);
  for (int i = 0; i < vec.size(); ++i)
    printf("%10.8g ", vec(i));
  printf("\n");
}

inline void PrintPretty(const char *identifier,
                        const Vector<std::complex<float>> &vec) {
  printf("%s:\n", identifier);
  for (int i = 0; i < vec.size(); ++i)
    printf("%10.8g%-+8.8gj ", vec(i).real(), vec(i).imag());
  printf("\n");
}

inline void PrintPretty(const char *identifier,
                        const Vector<std::complex<double>> &vec) {
  printf("%s:\n", identifier);
  for (int i = 0; i < vec.size(); ++i)
    printf("%10.8g%-+8.8gj ", vec(i).real(), vec(i).imag());
  printf("\n");
}

} // namespace lila
