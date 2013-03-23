//
//  nrm2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _nrm2_h
#define _nrm2_h

/// @file nrm2.h Computes the Euclidean norm of a vector, taking care not to cause unnecessary overflow.

#include <cmath>
#include "latl.h"

namespace LATL
{
   /// @brief Computes the Euclidean norm of a real vector.
   /// @return Euclidean norm ||x|| = sqrt(x'.x).
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of vector x.
   /// @ingroup VEC
   
   template <typename real_t>
   real_t NRM2(int_t n, real_t *x, int_t incx)
   {
      using std::abs;
      using std::sqrt;
      const real_t one(1.0);
      const real_t zero(0.0);
      real_t norm(0.0);
      real_t scale,ssq,a,b;
      int_t i;
      if((n>0)&&(incx>0))
      {
         if(n==1)
         {
            norm=abs(x[0]);         
         }
         else
         {
            scale=zero;
            ssq=one;
            for(i=0;i<n*incx;i+=incx)
            {
               if(x[i]!=zero)
               {
                  a=abs(x[i]);
                  if(scale<a)
                  {
                     b=scale/a;
                     ssq=one+ssq*(b*b);
                     scale=a;
                  }
                  else
                  {
                     b=a/scale;
                     ssq+=(b*b);
                  }
               }
            }
            norm=scale*sqrt(ssq);
         }
      }
      return norm;
   }
   
   /// @brief Computes the Euclidean norm of a complex vector.
   /// @return Euclidean norm ||x|| = sqrt(x'.x).
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of vector x.
   /// @ingroup VEC

   template <typename real_t>
   real_t NRM2(int_t n, complex<real_t> *x, int_t incx)
   {
      using std::abs;
      using std::real;
      using std::imag;
      using std::sqrt;
      const real_t one(1.0);
      const real_t zero(0.0);
      real_t norm(0.0);
      real_t scale,ssq,a,b;
      int_t i;
      if((n>0)&&(incx>0))
      {
         scale=zero;
         ssq=one;
         for(i=0;i<n*incx;i+=incx)
         {
            if(real(x[i])!=zero)
            {
               a=abs(real(x[i]));
               if(scale<a)
               {
                  b=scale/a;
                  ssq=one+ssq*(b*b);
                  scale=a;
               }
               else
               {
                  b=a/scale;
                  ssq+=(b*b);
               }
            }
            if(imag(x[i])!=zero)
            {
               a=abs(imag(x[i]));
               if(scale<a)
               {
                  b=scale/a;
                  ssq=one+ssq*(b*b);
                  scale=a;
               }
               else
               {
                  b=a/scale;
                  ssq+=(b*b);
               }
            }
         }
         norm=scale*sqrt(ssq);
      }
      return norm;
   }
   
}
#endif
