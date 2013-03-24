//
//  lassq.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 6/1/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lassq_h
#define _lassq_h

/// @file lassq.h Calculates the sum of the squares of elements of a vector.

#include <cmath>
#include "latl.h"

namespace LATL
{
   
   /// @brief Calculates the sum of the squares of elements of a real-valued vector.
   ///
   /// LASSQ() returns the values scl and ssq such that
   ///
   ///      (scl^2)*ssq = (x0)^2 + (x1)^2 + (x2)^2 + ... + (x[n-1])^2 + (scale)^2*sumsq
   ///
   /// where scale and sumsq are input values subsequently overwritten by scl and ssq.
   ///
   /// @tparam real_t Floating point type.
   /// @param n Number of elements of array X to be included.
   /// @param X Pointer to the first element to be included.
   /// @param incx Increment between successive elements of X.
   /// @param scale On entry, the value scale in the equation above.
   /// On exit, the scaling factor for the sum of squares.
   /// @param sumsq On entry, the value sumsq in the equation above.
   /// On exit, the basic sum of the squares from which scl has been factored out.  Assumed to be non-negative.
   /// @ingroup VEC
   
   template< typename real_t>
   void LASSQ( const int_t n, real_t * const X, const int_t incx, real_t &scale, real_t &sumsq)
   {
      using std::abs;
      using std::isnan;
      if (n > 0)
      {
         real_t zero(0.0), absxi;
         real_t temp=zero;
         for (int_t ix = 0; ix < n*incx; ix += incx)
         {
            absxi = abs(X[ix]);
            if (absxi > zero)
            {
               if (scale < absxi)
               {
                  temp = scale/absxi;
                  sumsq = 1 + sumsq*(temp*temp);
                  scale = absxi;
               }
               else
               {
                  temp = absxi/scale;
                  sumsq += (temp*temp);
               }
            }
            else if (isnan(absxi))
            {
               sumsq = temp;
            }
         }
      }
      return;
   }
   
   
   /// @brief Calculates the sum of the squares of elements of a complex-valued vector.
   ///
   /// LASSQ() returns the values scl and ssq such that
   ///
   ///      (scl^2)*ssq = (x0)^2 + (x1)^2 + (x2)^2 + ... + (x[n-1])^2 + (scale)^2*sumsq
   ///
   /// where scale and sumsq are input values subsequently overwritten by scl and ssq.
   ///
   /// @tparam real_t Floating point type.
   /// @param n Number of elements of array X to be included.
   /// @param X Pointer to the first (complex) element to be included.
   /// @param incx Increment between successive elements of X.
   /// @param scale On entry, the value scale in the equation above, assumed to be non-negative.
   /// On exit, the scaling factor for the sum of squares = maximum value among scale, abs(real(x[i])), and abs(imag(x[i])).
   /// @param sumsq On entry, the value sumsq in the equation above, assumed to be at least unity.
   /// On exit, the value ssq, which satisfies 1 <= ssq <= (sumsq + 2*n).
   /// @ingroup VEC
   
   template< typename real_t>
   void LASSQ( const int_t n, complex<real_t> * const X, const int_t incx, real_t &scale, real_t &sumsq)
   {
      using std::abs;
      using std::isnan;
      if (n > 0)
      {
         real_t zero(0.0);
         real_t temp=zero;
         real_t temp2=zero;
         for (int_t ix = 0; ix < n*incx; ix += incx)
         {
            temp = abs(real(X[ix]));
            if (temp>zero)
            {
               if (scale < temp)
               {
                  temp2 = scale/temp;
                  sumsq = 1 + sumsq*(temp2*temp2);
                  scale = temp;
               }
               else
               {
                  temp2 = temp/scale;
                  sumsq += (temp2*temp2);
               }
            }
            else if (isnan(temp))
            {
               sumsq = temp;
            }
            temp = abs(imag(X[ix]));
            if (temp>zero)
            {
               if (scale < temp)
               {
                  temp2 = scale/temp;
                  sumsq = 1+ sumsq*(temp2*temp2);
                  scale = temp;
               }
               else
               {
                  temp2 = temp/scale;
                  sumsq += temp2*temp2;
               }
            }
            else if (isnan(temp))
            {
               sumsq = temp;
            }
         }
      }
      return;
   }
}

#endif
