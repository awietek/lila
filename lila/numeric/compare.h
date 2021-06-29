#pragma once

#include <cmath>  
#include <cassert>

#include <lila/matrix.h>
#include <lila/vector.h>
#include <lila/numeric/complex.h>
#include <lila/numeric/precision.h>

namespace lila {
  
  template <class coeff_t>
  inline bool close(coeff_t const& x, coeff_t const& y,
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  { return (std::abs(x - y) <= ( atol + rtol * std::abs(y))); }

  
  template <class coeff_t>
  inline bool close(Vector<coeff_t> const& vec1,
		    Vector<coeff_t> const& vec2, 
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  {
    assert(vec1.size() == vec2.size());
    return std::equal(vec1.begin(), vec1.end(), vec2.begin(),
		      [&atol, &rtol](coeff_t x, coeff_t y) {
		      	return close(x, y, atol, rtol); 
		      });
  }
  
  template <class coeff_t>
  inline bool close(Matrix<coeff_t> const& mat1,
		    Matrix<coeff_t> const& mat2, 
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  {
    assert(mat1.nrows() == mat2.nrows());
    assert(mat1.ncols() == mat2.ncols());
    return std::equal(mat1.begin(), mat1.end(), mat2.begin(),
		      [&atol, &rtol](coeff_t x, coeff_t y) {
		      	return close(x, y, atol, rtol); 
		      });
  }

  template <class coeff_t>
  inline bool close(Vector<coeff_t> const& vec, coeff_t const& val, 
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  {
    return std::all_of(vec.begin(), vec.end(),
		       [&atol, &rtol, &val](coeff_t x) {
			 return close(x, val, atol, rtol); 
		       });
  }
  
  template <class coeff_t>
  inline bool close(Matrix<coeff_t> const& mat, coeff_t const& val, 
		    real_t<coeff_t> atol = lila::atol<coeff_t>::val(), 
		    real_t<coeff_t> rtol = lila::rtol<coeff_t>::val())
  {
    return std::all_of(mat.begin(), mat.end(), 
		       [&atol, &rtol, &val](coeff_t x) {
		      	return close(x, val, atol, rtol); 
		      });
  }
  

  template <class coeff_t, template<class> class object_t>
  inline bool equal(object_t<coeff_t> const& mat1,
		    object_t<coeff_t> const& mat2)
  {
    assert(mat1.size() == mat2.size());
    return std::equal(mat1.begin(), mat1.end(), mat2.begin());
  }
  
}
