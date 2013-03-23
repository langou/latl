//
//  sum1.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/20/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _sum1_h
#define _sum1_h

/// @file sum1.h Finds the sum of absolute values of complex vector elements.

#include <cmath>
#include "latl.h"

namespace LATL
{
   /// @brief Find the sum of the absolute values of complex vector elements.
   /// @return The sum of |x[i]| for i=0,...,n-1.
   /// @tparam real_t Floating point type.
   /// @param n Length of the vector x.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of the complex vector x; incx > 0.
   /// @ingroup VEC

   template <typename real_t>
   real_t SUM1(int_t n, complex<real_t> *x, int_t incx)
   {
      using std::abs;
      int_t i;
      real_t sum(0.0);
      if(n>0)
      {
         if(incx==1)
         {
            for(i=0;i<n;i++)
               sum+=abs(x[i]);
         }
         else if(incx>1)
         {
            for(i=0;i<n*incx;i+=incx)
               sum+=abs(x[i]);
         }
      }
      return sum;
   }
}
#endif
