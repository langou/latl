//
//  lapeq.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _lapeq_h
#define _lapeq_h

/// @file lapeq.h Performs the matrix operation A+=alpha*B.

#include "latl.h"

namespace LATL
{
   /// @brief Performs the matrix operation A+=alpha*B for real matrices A and B.
   /// @tparam real_t Floating point type.
   /// @param m The number of rows of the matrices A and B.  m>=0
   /// @param n The number of columns of the matrices A and B.  n>=0
   /// @param A Real m-by-n matrix.
   /// @param ldA Column length of the matrix A.  ldA>=m.
   /// @param alpha Real scalar.
   /// @param B Real m-by-n matrix.
   /// @param ldB Column length of the matrix B.  ldB>=m.

   template<typename real_t>
   void LAPEQ(int_t m, int_t n, real_t *A, int_t ldA,real_t alpha, real_t *B, int_t ldB)
   {
      const real_t one(1.0);
      if(alpha==one)
      {
         for(int_t j=0;j<n;j++)
         {
            for(int_t i=0;i<m;i++)
               A[i]+=B[i];
            A+=ldA;
            B+=ldB;
         }

      }
      else if(alpha==-one)
      {
         for(int_t j=0;j<n;j++)
         {
            for(int_t i=0;i<m;i++)
               A[i]-=B[i];
            A+=ldA;
            B+=ldB;
         }
      }
      else
      {
         for(int_t j=0;j<n;j++)
         {
            for(int_t i=0;i<m;i++)
               A[i]+=alpha*B[i];
            A+=ldA;
            B+=ldB;
         }
      }
   }

   /// @brief Performs the matrix operation A+=alpha*B for complex matrices A and B.
   /// @tparam real_t Floating point type.
   /// @param m The number of rows of the matrices A and B.  m>=0
   /// @param n The number of columns of the matrices A and B.  n>=0
   /// @param A Complex m-by-n matrix.
   /// @param ldA Column length of the matrix A.  ldA>=m.
   /// @param alpha Complex scalar.
   /// @param B Complex m-by-n matrix.
   /// @param ldB Column length of the matrix B.  ldB>=m.

   template<typename real_t>
   inline void LAPEQ(int_t m, int_t n, complex<real_t> *A, int_t ldA,complex<real_t> alpha, complex<real_t> *B, int_t ldB)
   {
      LAPEQ< complex<real_t> >(m,n,A,ldA,alpha,B,ldB);
   }
}

#endif
