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

#ifndef LILA_ALGORITHM_GRAMSCHMIDT_H_
#define LILA_ALGORITHM_GRAMSCHMIDT_H_

#include <vector>

namespace lila {
  template <class vector_t>
  std::vector<vector_t> gramschmidt(std::vector<vector_t>& ws)
  {
    auto es = ws;
    for (int i=0; i<(int)ws.size(); ++i)
      {
	auto v = ws[i];
	for (int k=0; k<i; ++k)
	  v -= Dot(es[k], ws[i])*es[k];
	v /= (typename vector_t::coeff_type)Norm(v);
	es[i] = v;
      }
    return es;   
  }

  template <class vector_t>
  bool has_full_rank(std::vector<vector_t>& ws)
  {
    int size = (int)ws.size();
    auto ovlps = lila::Zeros<typename vector_t::coeff_type>(size, size);
    for (int i=0; i<size; ++i)
      {
    	ovlps(i, i) = lila::Dot(ws[i], ws[i]);
    	for (int j=i+1; j<size; ++j)
    	  {
    	    ovlps(i, j) = lila::Dot(ws[i], ws[j]);
    	    ovlps(j, i) = lila::conj(ovlps(i, j));
    	  }
      }
    auto eigs = EigenvaluesH(ovlps);
    for (auto e : eigs)
      if (close(e, 0.)) return false;
    return true;
  }


}

#endif
