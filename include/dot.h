//
//  dot.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _dot_h
#define _dot_h

/// @file dot.h Computes dot products of two vectors.

#include "latl.h"

namespace LATL
{
   
   /// @brief Forms the dot product of two real vectors.
   /// @return Dot product of real vectors x and y.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors x and y.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to real vector y.
   /// @param incy Increment of y.
   /// @ingroup BLAS

   template <typename real_t>
   real_t DOT(int_t n, real_t *x, int_t incx, real_t *y, int_t incy)
   {
      real_t sum(0.0);
      int_t i,ix,iy;
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
               sum+=x[i]*y[i];
         }
         else
         {
            ix = incx>0 ? 0 : (1-n)*incx;
            iy = incy>0 ? 0 : (1-n)*incy;
            for(i=0;i<n;i++)
            {
               sum+=x[ix]*y[iy];
               ix+=incx;
               iy+=incy;
            }
         }
      }
      return sum;
   }
      
   /// @brief Forms the dot product of the conjugate of a complex vector with a complex vector.
   ///
   /// Computes
   ///
   ///     conj(x).y
   /// for complex vectors x and y.
   /// @return Dot product of complex vectors conj(x) and y.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors x and y.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to complex vector y.
   /// @param incy Increment of y.
   /// @ingroup BLAS

   template <typename real_t>
   complex<real_t> DOTC(int_t n, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy)
   {
      using std::conj;
      complex<real_t> sum(0.0,0.0);
      int_t i,ix,iy;
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
               sum+=conj(x[i])*y[i];
         }
         else
         {
            ix = incx>0 ? 0 : (1-n)*incx;
            iy = incy>0 ? 0 : (1-n)*incy;
            for(i=0;i<n;i++)
            {
               sum+=conj(x[ix])*y[iy];
               ix+=incx;
               iy+=incy;
            }
         }
      }
      return sum;
   }
   
   /// @brief Forms the dot product of two complex vectors.
   /// @return Dot product of complex vectors x and y.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors x and y.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to complex vector y.
   /// @param incy Increment of y.
   /// @ingroup BLAS

   template <typename real_t>
   complex<real_t> DOT(int_t n, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy)
   {
      complex<real_t> sum(0.0,0.0);
      int_t i,ix,iy;
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
               sum+=x[i]*y[i];
         }
         else
         {
            ix = incx>0 ? 0 : (1-n)*incx;
            iy = incy>0 ? 0 : (1-n)*incy;
            for(i=0;i<n;i++)
            {
               sum+=x[ix]*y[iy];
               ix+=incx;
               iy+=incy;
            }
         }
      }
      return sum;
   }
}
#endif
