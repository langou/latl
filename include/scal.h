//
//  scal.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/6/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _scal_h
#define _scal_h

/// @file scal.h Scales a vector by a scalar.

#include "latl.h"

namespace latl
{
   /// @brief Scales a real vector x by a real scalar alpha.
   ///
   ///        x := alpha*x
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x.
   /// @param alpha Real scalar.
   /// @param x Pointer to real vector.
   /// @param incx Increment of vector x.
   /// @ingroup VEC

   template <typename real_t>
   void scal(int_t n, real_t alpha, real_t *x, int_t incx)
   {
      int_t i,ix;
      if((n>0)&&(incx>0))
      {
         if(incx==1)
         {
            for(i=0;i<n;i++)
               x[i]*=alpha;
         }
         else
         {
            ix=0;
            for(i=0;i<n;i++)
            {
               x[ix]*=alpha;
               ix+=incx;
            }
         }
      }
   }
   
   /// @brief Scales a complex vector x by a real scalar alpha.
   ///
   ///        x := alpha*x
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x.
   /// @param alpha Real scalar.
   /// @param x Pointer to complex vector.
   /// @param incx Increment of vector x.
   /// @ingroup VEC

   template <typename real_t>
   void scal(int_t n, real_t alpha, complex<real_t> *x, int_t incx)
   {
      int_t i,ix;
      if((n>0)&&(incx>0))
      {
         if(incx==1)
         {
            for(i=0;i<n;i++)
               x[i]*=alpha;
         }
         else
         {
            ix=0;
            for(i=0;i<n;i++)
            {
               x[ix]*=alpha;
               ix+=incx;
            }
         }
      }
   }

   /// @brief Scales a complex vector x by a complex scalar alpha.
   ///
   ///        x := alpha*x
   /// @tparam real_t Floating point type.
   /// @param n Length of vector x.
   /// @param alpha Complex scalar.
   /// @param x Pointer to complex vector.
   /// @param incx Increment of vector x.
   /// @ingroup VEC

   template <typename real_t>
   void scal(int_t n, complex<real_t> alpha, complex<real_t> *x, int_t incx)
   {
      int_t i,ix;
      if((n>0)&&(incx>0))
      {
         if(incx==1)
         {
            for(i=0;i<n;i++)
               x[i]*=alpha;
         }
         else
         {
            ix=0;
            for(i=0;i<n;i++)
            {
               x[ix]*=alpha;
               ix+=incx;
            }
         }
      }
   }
}
#endif
