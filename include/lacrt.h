//
//  lacrt.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/24/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lacrt_h
#define _lacrt_h

/// @file lacrt.h Modified rotation of complex pairs of points.

#include "latl.h"

namespace LATL
{
   /// @brief Computes a modified rotation of complex pairs of points.
   /// 
   /// For complex scalars c and s and complex vectors x and y
   ///
   ///        x := cx+sy
   ///        y := cy-sx
   /// is computed pointwise for each element of x and y.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors x and y.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to complex vector y.
   /// @param incy Increment of vector y.
   /// @param c Complex scalar.
   /// @param s Complex scalar.
   /// @ingroup AUX
   
   template<typename real_t>
   void LACRT(int_t n, complex<real_t> *x, int_t incx, complex<real_t> *y,int_t incy,complex<real_t> c, complex<real_t> s)
   {
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(int_t i=0;i<n;i++)
            {
               complex<real_t> temp=c*x[i]+s*y[i];
               y[i]=c*y[i]-s*x[i];
               x[i]=temp;
            }
         }
         else
         {
            int_t ix=(incx>0)?0:(1-n)*incx;
            int_t iy=(incy>0)?0:(1-n)*incy;
            for(int_t i=0;i<n;i++)
            {
               complex<real_t> temp=c*x[ix]+s*y[iy];
               y[iy]=c*y[iy]-s*x[ix];
               x[ix]=temp;
               ix+=incx;
               iy+=incy;
            }            
         }
      }
   }
}
#endif
