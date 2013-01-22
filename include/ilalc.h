//
//  ilalc.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 5/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _ilalc_h
#define _ilalc_h

/// @file ilalc.h Scans matrix for its last nonzero column.

#include "latl.h"

namespace latl
{
   /// @brief Scans real matrix for its last nonzero column.
   /// @return The index of the last nonzero column in [0,n-1].  If A is the zero matrix, -1 is returned.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows in A; m >= 0.
   /// @param n Number of columns in A; n >= 0.
   /// @param A Pointer to real matrix A. 
   /// @param ldA Column length of A; ldA >= m.
   /// @ingroup MAT

   template<typename real_t>
   int_t ilalc(const int_t m,const int_t n,real_t *A,const int_t ldA)
   {
      const real_t zero(0.0);
      if(n>0)
      {
         A+=n*ldA;
         for(int_t j=n-1;j>=0;--j)
         {
            A-=ldA;
            for(int_t i=0;i<m;i++)
               if(A[i]!=zero)
                  return j;
         }
      }
      return -1;
   }

   /// @brief Scans complex matrix for its last nonzero column.
   /// @return The index of the last nonzero column in [0,n-1].  If A is the zero matrix, -1 is returned.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows in A; m >= 0.
   /// @param n Number of columns in A; n >= 0.
   /// @param A Pointer to complex matrix A. 
   /// @param ldA Column length of A; ldA >= m.
   /// @ingroup MAT

   template<typename real_t>
   int_t ilalc(const int_t m,const int_t n,complex<real_t> *A,const int_t ldA)
   {
      const complex<real_t> zero(0.0,0.0);
      if(n>0)
      {
         A+=n*ldA;
         for(int_t j=n-1;j>=0;--j)
         {
            A-=ldA;
            for(int_t i=0;i<m;i++)
               if(A[i]!=zero)
                  return j;
         }
      }
      return -1;
   }
}

#endif
