//
//  asum.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _asum_h
#define _asum_h

/// @file asum.h Finds the sum of absolute values of vector elements.

#include <cmath>
#include "latl.h"

namespace LATL
{
   /// @brief Find the sum of the absolute values of real vector elements.
   /// @return The sum of |x[i]| for i=0,...,n-1.
   /// @tparam real_t Floating point type.
   /// @param n Length of the vector x.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of the real vector x; incx > 0.
   /// @ingroup BLAS
   
   template <typename real_t>
   real_t ASUM(int_t n, real_t *x, int_t incx)
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

   /// @brief Find the sum of the absolute values of the real and imaginary parts of complex vector elements.
   /// @return The sum of |real(x[i])|+|imag(x[i])| for i=0,...,n-1.
   /// @tparam real_t Floating point type.
   /// @param n Length of the vector x.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of the complex vector x; incx > 0.
   /// @ingroup BLAS
   
   template <typename real_t>
   real_t ASUM(int_t n, complex<real_t> *x, int_t incx)
   {
      using std::abs;
      using std::real;
      using std::imag;
      int_t i;
      real_t sum(0.0);
      if(n>0)
      {
         if(incx==1)
         {
            for(i=0;i<n;i++)
               sum+=abs(real(x[i]))+abs(imag(x[i]));
         }
         else if(incx>1)
         {
            for(i=0;i<n*incx;i+=incx)
               sum+=abs(real(x[i]))+abs(imag(x[i]));
         }
      }
      return sum;
   }
}
#endif
