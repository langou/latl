//
//  swap.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _swap_h
#define _swap_h

/// @file swap.h Interchanges two vectors.

#include "latl.h"

namespace latl
{
   /// @brief Interchanges two real vectors.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to real vector y.
   /// @param incy Incremant of vector y.
   /// @ingroup VEC

   template <typename real_t>
   void swap(int_t n, real_t *x, int_t incx, real_t *y, int_t incy)
   {
      real_t temp;
      int_t i,ix,iy;
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
            {
               temp=x[i];
               x[i]=y[i];
               y[i]=temp;
            }
         }
         else
         {
            ix = incx>0 ? 0 : (1-n)*incx;
            iy = incy>0 ? 0 : (1-n)*incy;
            for(i=0;i<n;i++)
            {
               temp=x[ix];
               x[ix]=y[iy];
               y[iy]=temp;
               ix+=incx;
               iy+=incy;
            }
         }
      }
   }
   
   /// @brief Interchanges two complex vectors.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to complex vector y.
   /// @param incy Incremant of vector y.
   /// @ingroup VEC

   template <typename real_t>
   void swap(int_t n, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy)
   {
      swap< complex<real_t> >(n,x,incx,y,incy);
   }
}

#endif
