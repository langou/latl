//
//  lacgv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/24/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lacgv_h
#define _lacgv_h

/// @file lacgv.h Conjugates a complex vector.

#include "latl.h"

namespace LATL
{
   /// @brief Conjugates a complex vector x of length n.
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of vector x.
   /// @ingroup BLAS

   template<typename real_t>
   void LACGV(int_t n,complex<real_t> *x,int_t incx)
   {
      using std::conj;

      if(incx==1)
      {
         for(int_t i=0;i<n;i++)
            x[i]=conj(x[i]);
      }
      else
      {
         int_t ix=(incx>0)?0:(1-n)*incx;
         for(int_t i=0;i<n;i++)
         {
            x[ix]=conj(x[ix]);
            ix+=incx;
         }
      }
   }
}

#endif
