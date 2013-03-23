//
//  rscl.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/21/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _rscl_h
#define _rscl_h

/// @file rscl.h Multiplies a vector by the reciprocal of a scalar.

#include <limits>
#include "labad.h"
#include "scal.h"
#include "latl.h"

namespace LATL
{
   /// @brief Multiplies a real vector x by the real scalar 1/a.
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x.
   /// @param a Real scalar; vector x is multiplied by 1/a.
   /// @param x Pointer to real vector of length n.
   /// @param incx Increment of vector x.
   /// @ingroup VEC
   
   template<typename real_t> 
   void rscl(int_t n, real_t a, real_t *x, int_t incx)
   {
      using std::numeric_limits;
      using LATL::LABAD;
      using LATL::SCAL;
      const real_t one(1.0);
      const real_t zero(0.0);
      real_t a_inv;
      
      if(n>0)
      {
         real_t small_num=numeric_limits<real_t>::min();
         real_t big_num=one/small_num;
         LABAD<real_t>(small_num,big_num);
         real_t den=a;
         real_t num=one;
         bool done=0;
         while(!done)
         {
            real_t den1=den*small_num;
            real_t num1=num/big_num;
            if((abs(den1)>abs(num))&&(num!=zero))
            {
               a_inv=small_num;
               den=den1;
               done=0;
            }
            else if(abs(num1)>abs(den))
            {
               a_inv=big_num;
               num=num1;
               done=0;
            }
            else
            {
               a_inv=num/den;
               done=1;
            }
            SCAL(n,a_inv,x,incx);
         }
      }
   }

   /// @brief Multiplies a complex vector x by the real scalar 1/a.
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x.
   /// @param a Real scalar; vector x is multiplied by 1/a.
   /// @param x Pointer to complex vector of length n.
   /// @param incx Increment of vector x.
   /// @ingroup VEC
   
   template<typename real_t> 
   void rscl(int_t n, real_t a, complex<real_t> *x, int_t incx)
   {
      using std::numeric_limits;
      using LATL::LABAD;
      using LATL::SCAL;
      const real_t one(1.0);
      const real_t zero(0.0);
      real_t a_inv;
      
      if(n>0)
      {
         real_t small_num=numeric_limits<real_t>::min();
         real_t big_num=one/small_num;
         LABAD<real_t>(small_num,big_num);
         real_t den=a;
         real_t num=one;
         bool done=0;
         while(!done)
         {
            real_t den1=den*small_num;
            real_t num1=num/big_num;
            if((abs(den1)>abs(num))&&(num!=zero))
            {
               a_inv=small_num;
               den=den1;
               done=0;
            }
            else if(abs(num1)>abs(den))
            {
               a_inv=big_num;
               num=num1;
               done=0;
            }
            else
            {
               a_inv=num/den;
               done=1;
            }
            SCAL(n,a_inv,x,incx);
         }
      }
   }
}

#endif
