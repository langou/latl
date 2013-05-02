//
//  gehd2.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 2/27/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _gehd2_h
#define _gehd2_h

/// @file gehd2.h Reduces a general matrix to upper Hessenberg form.

#include <cstddef>
#include <algorithm>
#include "larfg.h"
#include "larf.h"
#include "latl.h"

namespace LATL
{
   /// @brief Reduces a general real matrix to upper Hessenberg form.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param A Real n-by-n matrix.  On exit, the upper triangle and the first
   /// subdiagonal of A are overwritten with the upper Hessenberg matrix H, and
   /// the elements below the first subdiagonal, with the array tau, represent
   /// the orthogonal matrix Q as a product of elementary reflectors.
   /// @param ldA Column length of matrix A.  ldA>=n
   /// @param tau Vector of length n-1.  The scalar factors of the elementary reflector.
   /// @ingroup COMP

   template <typename real_t>
   int GEHD2(int_t n, real_t *A, int_t ldA, real_t *tau)
   {
      const real_t one(1.0);
      using std::min;
      if(n<0)
         return -1;
      else if(ldA<n)
         return -3;

      real_t *w = new real_t[n];
      real_t *B=A+ldA;
      for(int_t i=0;i<n-1;i++)
      {
         real_t alpha=A[i+1];
         real_t *v=A+min(i+2,n-1);
         LARFG(n-i-1,alpha,v,1,tau[i]);
         A[i+1]=one;
         LARF('R',n,n-i-1,A+i+1,1,tau[i],B,ldA,w);
         LARF('L',n-i-1,n-i-1,A+i+1,1,tau[i],B+i+1,ldA,w);
         A[i+1]=alpha;
         A+=ldA;
         B+=ldA;
      }
      delete [] w;
      return 0;
   }
}
#endif
