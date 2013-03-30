//
//  dotx.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _dotx_h
#define _dotx_h

/// @file dotx.h Computes dot products of two vectors in mixed precision.

#include "latl.h"

namespace LATL
{
   /// @brief Forms the dot product of two real vectors using extended precision for accumulation.
   ///
   /// Computes
   ///
   ///     b+x.y 
   /// using extended precision for the dot product accumulation; b is a scalar, x and y are vectors.
   /// @return Dot product of real vectors x and y plus b.
   /// @tparam real_t Floating point type.
   /// @tparam xreal_t Real floating point type of higher precision.
   /// @param n Length of vectors x and y.
   /// @param b Real scalar, initializes dot product accumulation.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to real vector y.
   /// @param incy Increment of y.
   /// @ingroup BLAS

   template <typename real_t,typename xreal_t>
   real_t DOTX(int_t n, real_t b, real_t *x, int_t incx, real_t *y, int_t incy)
   {
      xreal_t sum(b);
      real_t sumf;
      int_t i,ix,iy;
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
               sum+=((xreal_t)x[i])*((xreal_t)y[i]);
         }
         else
         {
            ix = (incx>0)?0:(1-n)*(incx);
            iy = (incy>0)?0:(1-n)*(incy);
            for(i=0;i<n;i++)
            {
               sum+=((xreal_t)x[ix])*((xreal_t)y[iy]);
               ix+=incx;
               iy+=incy;
            }
         }
      }
      sumf=sum;
      return sumf;
   }

   
   /// @brief Forms the dot product of two real vectors using extended precision for accumulation.
   /// @return Dot product of real vectors x and y in extended precision.
   /// @tparam real_t Floating point type.
   /// @tparam xreal_t Real floating point type of higher precision.
   /// @param n Length of vectors x and y.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to real vector y.
   /// @param incy Increment of y.
   /// @ingroup BLAS

   template <typename real_t, typename xreal_t>
   xreal_t DOTX(int_t n, real_t *x, int_t incx, real_t *y, int_t incy)
   {
      xreal_t sum(0.0);
      int_t i,ix,iy;
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
               sum+=((xreal_t)x[i])*((xreal_t)y[i]);
         }
         else
         {
            ix = (incx>0)?0:(1-n)*(incx);
            iy = (incy>0)?0:(1-n)*(incy);
            for(i=0;i<n;i++)
            {
               sum+=((xreal_t)x[ix])*((xreal_t)y[iy]);
               ix+=incx;
               iy+=incy;
            }
         }
      }
      return sum;
   }
}
#endif
