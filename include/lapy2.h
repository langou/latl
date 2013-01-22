//
//  lapy2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 5/13/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lapy2_h
#define _lapy2_h

/// @file lapy2.h Finds sqrt(x^2+y^2), taking care not to cause unnecessary overflow.

#include <cmath>

namespace latl
{
   
   /// @brief Finds sqrt(x^2+y^2), taking care not to cause unnecessary overflow.
   /// @return sqrt(x^2+y^2)
   /// @tparam real_t Floating point type.
   /// @param x scalar value x
   /// @param y scalar value y
   /// @ingroup SCAL
   
   template<typename real_t>
   real_t lapy2(real_t x, real_t y)
   {
      using std::abs;
      using std::sqrt;
      
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t w,xabs,yabs,z,r;

      xabs=abs(x);
      yabs=abs(y);
      if(xabs>yabs)
      {
         w=xabs;
         z=yabs;
      }
      else
      {
         w=yabs;
         z=xabs;
      }
      if(z==zero)
         r=w;
      else
         r=w*sqrt(one+(z/w)*(z/w));
      return r;
   }
}

#endif
