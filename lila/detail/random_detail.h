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

#ifndef LILA_DETAIL_RANDOM_DETAIL_H_
#define LILA_DETAIL_RANDOM_DETAIL_H_

#include <random>
#include <complex>

namespace lila {
  namespace detail {

    template <class random_engine_t, class distribution_t>
    inline void random_number(random_engine_t& engine, 
			      distribution_t& distribution, float& number)
    { number = distribution(engine); }
    
    template <class random_engine_t, class distribution_t>
    inline void random_number(random_engine_t& engine, 
			      distribution_t& distribution, double& number)
    { number = distribution(engine); }

    template <class random_engine_t, class distribution_t>
    inline void random_number(random_engine_t& engine, 
			      distribution_t& distribution,
			      std::complex<float>& number)
    { number = {distribution(engine), distribution(engine)}; }

    template <class random_engine_t, class distribution_t>
    inline void random_number(random_engine_t& engine, 
			      distribution_t& distribution,
			      std::complex<double>& number)
    { number = {distribution(engine), distribution(engine)}; }
   
  }  
}

#endif  // LILA_DETAIL_RANDOM_DETAIL_H_
