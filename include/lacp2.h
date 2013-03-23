//
//  lacp2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lacp2_h
#define _lacp2_h

/// @file lacp2.h Copies a real matrix to a complex matrix.

#include "latl.h"

namespace LATL
{
   /// @brief Copies a real A to complex matrix B where A is either a full, upper triangular or lower triangular matrix.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the source matrix is upper triangular, lower triangular, or full:
   ///
   ///             if uplo = 'U' or 'u' then A is upper triangular
   ///             if uplo = 'L' or 'l' then A is lower triangular
   ///             otherwise, A is assumed to be full.
   ///
   /// @param m The number of rows of the matrix A. m>=0
   /// @param n The number of columns of the matrix A.  n>=0
   /// @param A Pointer to the real matrix A.
   /// @param ldA Column length of the matrix A.  ldA>=m
   /// @param B Pointer to the complex matrix B.
   /// @param ldB Column length of the matrix B.  ldB>=m
   /// @ingroup MAT
   
   template<typename real_t>
   void lacp2(char uplo,int_t m,int_t n,real_t *A,int_t ldA,complex<real_t> *B,int_t ldB)
   {
      switch(uplo)
      {
         case 'U':
         case 'u':
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;(i<m)&&(i<=j);i++)
                  B[i]=A[i];
               B+=ldB;
               A+=ldA;
            }
            break;
            
         case 'L':
         case 'l':
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=j;i<m;i++)
                  B[i]=A[i];
               B+=ldB;
               A+=ldA;
            }
            break;
            
         default:
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;i<m;i++)
                  B[i]=A[i];
               B+=ldB;
               A+=ldA;
            }
            break;
      }
   }
}

#endif
