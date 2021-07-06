#include "../catch.hpp"

#include <lila/all.h>
#include <unistd.h>

using namespace lila;

template <class coeff_t> MatrixView<coeff_t> return_matrix_view() {
  auto mat_del = Random<coeff_t>(5, 5);
  LilaPrint(mat_del);
  auto submat_del = mat_del({2, 4}, {2, 4});
  LilaPrint(mat_del.use_count());
  LilaPrint(submat_del.use_count());
  return submat_del;
}

template <class coeff_t> void test_views() {

  auto vec = Random<coeff_t>(8);
  auto subvec = vec({2,8,2});
  // LilaPrint(vec);
  // LilaPrint(subvec);

  // Check reference counting
  auto mat = Random<coeff_t>(5, 5);
  // LilaPrint(mat.use_count());
  {
    auto submat = mat({2, 4}, {2, 4});
    REQUIRE(submat.use_count() == 2);
    // LilaPrint(mat.use_count());
    // LilaPrint(submat.use_count());
    // LilaPrint(mat);
    // LilaPrint(submat);
  }
  REQUIRE(mat.use_count() == 1);

  // LilaPrint(mat.use_count());
  // printf("\n\n\n");


  // // check deletion of creating Matrix doesn't delete storage
  // auto view = return_matrix_view<coeff_t>();
  // LilaPrint(view);
  // LilaPrint(view.use_count());

  // printf("Setting submatrix to one\n");
  mat({2,4}, {2,4}) = (coeff_t)1.;
  REQUIRE(mat.use_count() == 1);
  // LilaPrint(mat);
  // LilaPrint(mat.use_count());
  // printf("\n\n\n");

  // printf("Adding two to submatrix \n");
  // mat({2,4}, {2,4}) += (coeff_t)2.;
  REQUIRE(mat.use_count() == 1);
  // LilaPrint(mat);
  // LilaPrint(mat.use_count());
  // printf("\n\n\n");

  // printf("subtracting five from submatrix \n");
  // mat({2,4}, {2,4}) -= (coeff_t)5.;
  // REQUIRE(mat.use_count() == 1);
  // LilaPrint(mat);
  // LilaPrint(mat.use_count());
  // printf("\n\n\n");



  // printf("setting first column equal to second column \n");
  // mat(ALL, 1) = mat(ALL, 0);
  // REQUIRE(mat.use_count() == 1);
  // LilaPrint(mat);
  // LilaPrint(mat.use_count());
  // printf("\n\n\n");

  // printf("setting first row to 0 \n");
  // mat(0, ALL) = Zeros<coeff_t>(5);
  // LilaPrint(mat);
  // LilaPrint(mat.use_count());
  // printf("\n\n\n");


  // // Check if storage isn't copied for views and 
  // // finally properly released
  // {
  //   auto mat = Zeros<coeff_t>(50000, 50000);
  //   printf("initial matrix allocated ....\n");
  //   sleep(10);
  //   {
  //     auto submat = mat({1, 9999}, {1, 9999});
  //     printf("submatrix view generated ....\n");
  //     sleep(10);
  //   }
  //   printf("submatrix view out of scope ....\n");
  //   sleep(10);
  // }
  // printf("initial matrix  out of scope ....\n");
  // sleep(10);


  // printf("---------------------------\n\n\n");

}

TEST_CASE("views", "[core]") {
  test_views<float>();
  test_views<double>();
  test_views<std::complex<float>>();
  test_views<std::complex<double>>();
}
