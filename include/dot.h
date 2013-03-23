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
   /// @ingroup VEC

   template <typename real_t>
   real_t dot(int_t n, real_t *x, int_t incx, real_t *y, int_t incy)
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
   /// @ingroup VEC

   template <typename real_t,typename xreal_t>
   real_t dotx(int_t n, real_t b, real_t *x, int_t incx, real_t *y, int_t incy)
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
   /// @ingroup VEC

   template <typename real_t, typename xreal_t>
   xreal_t dotx(int_t n, real_t *x, int_t incx, real_t *y, int_t incy)
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
   /// @ingroup VEC

   template <typename real_t>
   complex<real_t> dotc(int_t n, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy)
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
   /// @ingroup VEC

   template <typename real_t>
   complex<real_t> dot(int_t n, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy)
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
