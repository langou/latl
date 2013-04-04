//
//  lamch.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 4/4/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _lamch_h
#define _lamch_h

/// @file lamch.h Returns floating point machine parameters.

#include <limits>
#include <cctype>

namespace LATL
{
   /// @brief Returns floating point machine parameters.
   /// @tparam real_t Floating point type.
   /// @return Machine parameter specified by opt.
   /// @return Zero if input parameter not recognized.
   /// @param opt Specifies machine parameter to return:
   ///
   ///        'E' = relative machine precision
   ///        'P' = machine precision (epsilon)
   ///        'S' = safe minimum (sfmin), such that 1/sfmin does not overflow
   ///        'U' = underflow threshold
   ///        'O' = overflow threshold
   /// @ingroup AUX
   
   template <typename real_t>
   real_t LAMCH(char OPT)
   {
      using std::toupper;
      using std::numeric_limits;
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t prec=numeric_limits<real_t>::epsilon();
      const real_t rnd=numeric_limits<real_t>::round_error();
      const real_t eps=prec*rnd;
      const real_t maxval=numeric_limits<real_t>::max();
      const real_t minval=numeric_limits<real_t>::min();
      const real_t small=one/maxval;
      const real_t sfmin=(small>minval) ? small*(one+eps) : minval;
      const real_t underflow=numeric_limits<real_t>::min();
      const real_t overflow=numeric_limits<real_t>::max();

      const char opt=toupper(OPT);

      if(opt=='E')
         return eps;
      else if(opt=='S')
         return sfmin;
      else if(opt=='P')
         return prec;
      else if(opt=='U')
         return underflow;
      else if(opt=='O')
         return overflow;
      else
         return zero;
   }
}
#endif

