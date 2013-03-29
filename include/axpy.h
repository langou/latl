//
//  axpy.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/7/11.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _axpy_h
#define _axpy_h

/// @file axpy.h Computes a constant times a vector plus a vector.

#include "latl.h"

namespace LATL
{
   
   /// @brief Computes a real constant times a real vector plus a real vector.
   /// 
   /// For real vectors x and y and real scalar a, computes
   ///
   ///     y := ax + y
   ///
   /// @tparam real_t Floating point type.
   /// @param n Length of the vectors x and y.
   /// @param a Real constant.
   /// @param x Real vector of length n. 
   /// @param incx Increment of vector x.
   /// @param y Real vector of length n.  On exit, y contains the product ax+y.
   /// @param incy Increment of vector y.
   /// @ingroup BLAS

   template <typename real_t>
   void AXPY(int_t n, real_t a, real_t *x, int_t incx, real_t *y, int_t incy)
   {
      const real_t zero(0.0);
      int_t i,ix,iy;
      
      if((n>0)&&(a!=zero))
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
               y[i]+=a*x[i];
         }
         else
         {
            ix=(incx>0)?0:(1-n)*incx;
            iy=(incy>0)?0:(1-n)*incy;
            for(i=0;i<n;i++)
            {
               y[iy]+=a*x[ix];
               ix+=incx;
               iy+=incy;
            }
         }
      }
   }
   
   /// @brief Computes a complex constant times a complex vector plus a complex vector.
   /// 
   /// For complex vectors x and y and complex scalar a, computes
   ///
   ///     y := ax + y
   ///
   /// @tparam real_t Floating point type.
   /// @param n Length of the vectors x and y.
   /// @param a Complex constant.
   /// @param x Complex vector of length n. 
   /// @param incx Increment of vector x.
   /// @param y Complex vector of length n.  On exit, y contains the product ax+y.
   /// @param incy Increment of vector y.
   /// @ingroup BLAS

   template <typename real_t>
   void AXPY(int_t n, complex<real_t> a, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy)
   {
      const complex<real_t> zero(0.0,0.0);
      int_t i,ix,iy;
      
      if((n>0)&&(a!=zero))
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
               y[i]+=a*x[i];
         }
         else
         {
            ix=(incx>0)?0:(1-n)*incx;
            iy=(incy>0)?0:(1-n)*incy;
            for(i=0;i<n;i++)
            {
               y[iy]+=a*x[ix];
               ix+=incx;
               iy+=incy;
            }
         }
      }
   }
}
#endif
