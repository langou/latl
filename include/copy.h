//
//  copy.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _copy_h
#define _copy_h

/// @file copy.h Copies a vector to another vector.

#include "latl.h"

namespace LATL
{
   /// @brief Copies a real vector x to a real vector y, performing the operation y := x.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors x and y.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to real vector y.
   /// @param incy Increment of vector y.
   /// @ingroup BLAS

   template <typename real_t>
   void COPY(int_t n, real_t *x, int_t incx, real_t *y, int_t incy)
   {
      int_t i,ix,iy;
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
               y[i]=x[i];
         }
         else
         {
            ix = incx>0 ? 0 : (1-n)*incx;
            iy = incy>0 ? 0 : (1-n)*incy;
            for(i=0;i<n;i++)
            {
               y[iy]=x[ix];
               ix+=incx;
               iy+=incy;
            }
         }
      }
   }
   
   /// @brief Copies a complex vector x to a complex vector y, performing the operation y := x.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors x and y.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to complex vector y.
   /// @param incy Increment of vector y.
   /// @ingroup BLAS

   template <typename real_t>
   void COPY(int_t n, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy)
   {
      COPY< complex<real_t> >(n,x,incx,y,incy);
   }
}
#endif
