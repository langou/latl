//
//  lapy3.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/10/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lapy3_h
#define _lapy3_h

/// @file lapy3.h Finds sqrt(x^2+y^2+z^2), taking care not to cause unnecessary overflow.

#include <cmath>
#include <algorithm>

namespace latl
{
   /// @brief Finds sqrt(x^2+y^2+z^2), taking care not to cause unnecessary overflow.
   /// @return sqrt(x^2+y^2+z^2)
   /// @tparam real_t Floating point type.
   /// @param x scalar value x
   /// @param y scalar value y
   /// @param z scalar value z
   /// @ingroup SCAL

   template<typename real_t>
   real_t lapy3(real_t x, real_t y, real_t z)
   {
      using std::abs;
      using std::sqrt;
      using std::max;
      
      const real_t zero(0.0);
      real_t w,xabs,yabs,zabs,ret;
      xabs=abs(x);
      yabs=abs(y);
      zabs=abs(z);
      w=max(xabs,max(yabs,zabs));
      if(w==zero)
      {
         ret=xabs+yabs+zabs;
      }
      else
      {
         xabs=xabs/w;
         yabs=yabs/w;
         zabs=zabs/w;
         ret=w*sqrt(xabs*xabs+yabs*yabs+zabs*zabs);
      }
      return ret;
   }
}
#endif
