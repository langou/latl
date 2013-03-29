//
//  imax.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _imax_h
#define _imax_h

/// @file imax.h Finds the index of the maximal element of a vector.

#include <cmath>
#include "latl.h"

namespace LATL
{
   
   /// @brief Finds the index of the maximal element of a real vector.
   /// @return The index i such that |x[i]| is maximal, or -1 if either n<1 or incx<1.
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x; n>=0.
   /// @param x Pointer to real vector of length n;
   /// @param incx Increment of the real vector x; incx > 0
   /// @ingroup BLAS
   
   template <typename real_t>
   int_t IMAX(int_t n, real_t *x, int_t incx)
   {
      using std::abs;
      int_t i,ix,index;
      real_t maxval;
      if(n<1)
         index=-1;
      else if (incx<1)
         index=-1;
      else if(n==1)
         index=0;
      else if(incx==1)
      {
         maxval=abs(x[0]);
         index=0;
         for(i=1;i<n;i++)
         {
            if(abs(x[i])>maxval)
            {
               index=i;
               maxval=abs(x[i]);
            }
         }
      }
      else
      {
         maxval=abs(x[0]);
         index=0;
         ix=incx;
         for(i=1;i<n;i++)
         {
            if(abs(x[ix])>maxval)
            {
               index=i;
               maxval=abs(x[ix]);
            }
            ix+=incx;
         }
      }
      return index;
   }
   
   /// @brief Finds the index of the maximal element of a complex vector.
   /// @return The index i such that |real(x[i])|+|imag(x[i])| is maximal, or -1 if either n<1 or incx<1.
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x; n>=0.
   /// @param x Pointer to complex vector of length n;
   /// @param incx Increment of the real vector x; incx > 0
   /// @ingroup BLAS

   template <typename real_t>
   int_t IMAX(int_t n, complex<real_t> *x, int_t incx)
   {
      using std::abs;
      using std::real;
      using std::imag;
      
      int_t i,ix,index;
      real_t maxval,absval;
      if(n<1)
         index=-1;
      else if(incx<1)
         index=-1;
      else if(n==1)
         index=0;
      else if(incx==1)
      {
         absval=abs(real(x[0]))+abs(imag(x[0]));
         maxval=absval;
         index=0;
         for(i=1;i<n;i++)
         {
            absval=abs(real(x[i]))+abs(imag(x[i]));
            if(absval>maxval)
            {
               index=i;
               maxval=absval;
            }
         }
      }
      else
      {
         absval=abs(real(x[0]))+abs(imag(x[0]));
         maxval=absval;
         index=0;
         ix=incx;
         for(i=1;i<n;i++)
         {
            absval=abs(real(x[ix]))+abs(imag(x[ix]));
            if(absval>maxval)
            {
               index=i;
               maxval=absval;
            }
            ix+=incx;
         }
      }
      return index;
   }
}
#endif
