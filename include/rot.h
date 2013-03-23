//
//  rot.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _rot_h
#define _rot_h

/// @file rot.h Applies plane rotations to a series of points.

#include "latl.h"

namespace LATL
{
   /// @brief Applies real plane rotations to a series of real points (x,y).
   ///
   /// For a point (x,y) and real scalars c and s
   ///
   ///        [ x ] := [ c  s ] [ x ]
   ///        [ y ] := [-s  c ] [ y ]
   /// is computed.
   /// @tparam real_t Floating point type.
   /// @param n Number of points.
   /// @param x Pointer to vector of x coordinates.
   /// @param incx Increment of the vector x.
   /// @param y Pointer to vector of y coordinates.
   /// @param incy Increment of the vector y.
   /// @param c Real scalar.
   /// @param s Real scalar.
   /// @ingroup ROT
   
   template <typename real_t>
   void rot(int_t n, real_t *x, int_t incx, real_t *y, int_t incy, real_t c, real_t s)
   {
      real_t temp;
      int_t i,ix,iy;
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
            {
               temp=c*x[i]+s*y[i];
               y[i]=c*y[i]-s*x[i];
               x[i]=temp;
            }
         }
         else
         {
            ix = incx>0 ? 0 : (1-n)*incx;
            iy = incy>0 ? 0 : (1-n)*incy;
            for(i=0;i<n;i++)
            {
               temp=c*x[ix]+s*y[iy];
               y[iy]=c*y[iy]-s*x[ix];
               x[ix]=temp;
               ix+=incx;
               iy+=incy;
            }
         }
      }
   }
 
   /// @brief Applies complex plane rotations to a series of complex points (x,y).
   ///
   /// For a complex point (x,y), real scalar c and complex scalar s
   ///
   ///        [ x ] := [  c       s ] [ x ]
   ///        [ y ] := [ -conj(s) c ] [ y ]
   /// is computed.
   /// @tparam real_t Floating point type.
   /// @param n Number of points.
   /// @param x Pointer to vector of x coordinates.
   /// @param incx Increment of the vector x.
   /// @param y Pointer to vector of y coordinates.
   /// @param incy Increment of the vector y.
   /// @param c Real scalar.
   /// @param s Complex scalar.
   /// @ingroup ROT

   template <typename real_t>
   void rot(int_t n, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy, real_t c, complex<real_t> s)
   {
      using std::conj;
      complex<real_t> temp;
      int_t i,ix,iy;
      if(n>0)
      {
         if((incx==1)&&(incy==1))
         {
            for(i=0;i<n;i++)
            {
               temp=c*x[i]+s*y[i];
               y[i]=c*y[i]-conj(s)*x[i];
               x[i]=temp;
            }
         }
         else
         {
            ix = incx>0 ? 0 : (1-n)*incx;
            iy = incy>0 ? 0 : (1-n)*incy;
            for(i=0;i<n;i++)
            {
               temp=c*x[ix]+s*y[iy];
               y[iy]=c*y[iy]-conj(s)*x[ix];
               x[ix]=temp;
               ix+=incx;
               iy+=incy;
            }
         }
      }
   }
}
#endif
