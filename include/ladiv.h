//
//  ladiv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/12/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _ladiv_h
#define _ladiv_h

/// @file ladiv.h Performs complex division of two scalars taking care to avoid underflow.

#include <cmath>
#include "latl.h"

namespace latl
{
   /// @brief Performs complex division in real arithmetic.
   ///
   ///        p + iq = (a +ib) / (c + id)
   ///
   /// @tparam real_t Floating point type.
   /// @param a Real part of numerator.
   /// @param b Imaginary part of numerator.
   /// @param c Real part of denominator.
   /// @param d Imaginary part of denominator.
   /// @param p Real part of quotient.
   /// @param q Imaginary part of quotient.
   /// @ingroup SCAL

   template<typename real_t>
   void ladiv(real_t a, real_t b, real_t c, real_t d, real_t &p, real_t &q)
   {
      using std::abs;
      real_t e,f;
      if(abs(d)<abs(c))
      {
         e=d/c;
         f=c+d*e;
         p=(a+b*e)/f;
         q=(b-a*e)/f;
      }
      else
      {
         e=c/d;
         f=c+d*e;
         p=(b+a*e)/f;
         q=(-a+b*e)/f;
      }
   }
   
   /// @brief Performs complex division in real arithmetic with complex arguments.
   /// @return x/y
   /// @tparam real_t Floating point type.
   /// @param x Complex numerator.
   /// @param y Complex denominator.
   /// @ingroup SCAL

   template<typename real_t>
   complex<real_t> ladiv(complex<real_t> x, complex<real_t> y)
   {
      using std::abs;
      using std::real;
      using std::imag;
      real_t zr,zi;
      ladiv<real_t>(real(x),imag(x),real(y),imag(y),zr,zi);
      return complex<real_t>(zr,zi);
   }

}
#endif
