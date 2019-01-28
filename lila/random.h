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

#ifndef LILA_RANDOM_H_
#define LILA_RANDOM_H_

#include <random>

#include "complex.h"

#include "solve.h"
#include "detail/random_detail.h"

#include "matrix.h"
#include "vector.h"

namespace lila {
 
  template <class coeff_t>
  using uniform_dist_t = std::uniform_real_distribution<real_t<coeff_t>>; 

  template <class coeff_t>
  using normal_dist_t = std::normal_distribution<real_t<coeff_t>>; 

  template <class random_engine_t, class distribution_t, class coeff_t>
  class random_generator {
  public:
    random_generator(distribution_t& distribution, int seed)
      : distribution_(distribution)
    { engine_.seed(seed); }

    coeff_t operator() ()
    { 
      coeff_t number;
      detail::random_number(engine_, distribution_, number);
      return number; 
    }

    void seed(int seed) { engine_.seed(seed); }

  private:
    random_engine_t engine_;
    distribution_t distribution_;
  };  

  template <class coeff_t>
  using uniform_gen_t = random_generator<std::mt19937, uniform_dist_t<coeff_t>, 
					 coeff_t>; 
  template <class coeff_t>
  using normal_gen_t = random_generator<std::mt19937, normal_dist_t<coeff_t>, 
					 coeff_t>; 

  template <class matrix_t, class gen_t>
  void Random(matrix_t& mat, gen_t& gen, bool alter_generator=true)
  { 
    using coeff_t = typename matrix_t::coeff_type;
    if(alter_generator) std::for_each(mat.begin(), mat.end(), 
				      [&gen](coeff_t& c){ c = gen(); });
    else std::generate(mat.data(), mat.data() + mat.size(), gen); 
  }

  template <class coeff_t, class gen_t>
  Vector<coeff_t> Random(int m, gen_t& gen, bool alter_generator=true)
  { 
    Vector<coeff_t> vec(m);
    Random(vec, gen, alter_generator);
    return vec;
  }

  template <class coeff_t, class gen_t>
  Matrix<coeff_t> Random(int m, int n, gen_t& gen, bool alter_generator=true)
  { 
    Matrix<coeff_t> mat(m, n);
    Random(mat, gen, alter_generator);
    return mat;
  }

  template <class matrix_t, class gen_t>
  void RandomUnitary(matrix_t& mat, gen_t& gen, bool alter_generator=true)
  { 
    assert(mat.nrows() == mat.ncols()); // Could be dropped eventually
    using coeff_t = typename matrix_t::coeff_type;
    Random(mat, gen, alter_generator);
    std::vector<coeff_t> tau = QRDecompose(mat);
    mat = QRGetQ(mat, tau);    
    // TODO: get uniform Haar measure here !!!!
    // i.e. correct diagonal
    // cf. e.g. arxiv:math-ph/0609050v2
  }

}

#endif
