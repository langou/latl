//
//  ilalr.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 5/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _ilalr_h
#define _ilalr_h

/// @file ilalr.h Scans matrix for its last nonzero row.

#include <algorithm>
#include "latl.h"

namespace LATL
{
   /// @brief Scans real matrix for its last nonzero row.
   /// @return The index of the last nonzero row in [0,m-1].  If A is the zero matrix, -1 is returned.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows in A; m >= 0.
   /// @param n Number of columns in A; n >= 0.
   /// @param A Pointer to real matrix A. 
   /// @param ldA Column length of A; ldA >= m.
   /// @ingroup MAT

   template<typename real_t>
   int_t ilalr(int_t m, int_t n, real_t *A, int_t ldA)
   {
      using std::max;
      const real_t zero(0.0);
      int_t row=-1;
      if(m>0)
      {
         for(int_t j=0;j<n;j++)
         {
            int_t i=m-1;
            while((A[i]==zero)&&(i>=0))
               i--;
            row=max(row,i);
            A+=ldA;
         }
      }
      return row;
   }

   /// @brief Scans complex matrix for its last nonzero row.
   /// @return The index of the last nonzero row in [0,m-1].  If A is the zero matrix, -1 is returned.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows in A; m >= 0.
   /// @param n Number of columns in A; n >= 0.
   /// @param A Pointer to complex matrix A. 
   /// @param ldA Column length of A; ldA >= m.
   /// @ingroup MAT

   template<typename real_t>
   int_t ilalr(int_t m, int_t n, complex<real_t> *A, int_t ldA)
   {
      using std::max;
      const complex<real_t> zero(0.0,0.0);
      int_t row=-1;
      if(m>0)
      {
         for(int_t j=0;j<n;j++)
         {
            int_t i=m-1;
            while((A[i]==zero)&&(i>=0))
               i--;
            row=max(row,i);
            A+=ldA;
         }
      }
      return row;
   }
}

#endif
