//
//  labad.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/21/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _labad_h
#define _labad_h

/// @file labad.h Corrects underflow and overflow values.

#include <cmath>

namespace latl
{
   /// @brief Corrects underflow and overflow values.
   ///
   /// Takes as input the values computed by numeric_limits for underflow and
   /// overflow, and returns the square root of each of these values if the
   /// log of large is sufficiently big.  This function is intended to
   /// identify machines with a large exponent range and
   /// redefine the underflow and overflow limits to be the square roots of
   /// the values computed by numeric_limits.
   /// @tparam real_t Floating point type.
   /// @param[in,out] small Underflow threshold; if log10(large) is sufficiently
   /// big, sqrt(small) is returned, otherwise small is unchanged.
   /// @param[in,out] large Overflow threshold; if log10(large) is sufficiently
   /// big, sqrt(large) is returned, otherwise large is unchanged.
   /// @ingroup SCAL
   
   template<typename real_t>
   void labad(real_t &small, real_t &large)
   {
      using std::log10;
      using std::sqrt;
      const real_t big(2000.0);

      if(log10(large)>big)
      {
         small=sqrt(small);
         large=sqrt(large);
      }
   }
}

#endif
