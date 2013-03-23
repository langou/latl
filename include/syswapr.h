//
//  syswapr.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 7/9/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _syswapr_h
#define _syswapr_h

/// @file syswapr.h Applies an elementary permutation on the rows and columns of a symmetric matrix.

#include <cctype>
#include "swap.h"
#include "latl.h"

namespace LATL
{
   /// @brief Applies an elementary permutation on the rows and columns of a symmetric matrix.
   ///
   /// Given a pair of indices (i,j), with 0 <= i,j < n, the ith and jth rows and columns of the
   /// symmetric matrix A are interchanged.  A is stored as either upper or lower triangular, and
   /// only the upper or lower triangular part of the matrix is referenced.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the symmetric matrix A is stored as upper or lower triangular:
   ///
   ///        'U': A is stored as upper triangular; elements of A below the diagonal are not referenced.
   ///        'L': A is stored as lower triangular; elements of A above the diagonal are not referenced.
   /// @param n The order of the symmetric matrix A. n>=0
   /// @param A Pointer to real symmetric matrix [in/out].
   /// @param ldA Column length of matrix A. ldA>=n
   /// @param i First index for row and column transposition.  0<=i<n
   /// @param j Second index for row and column transposition. 0<=j<n
   
   template<typename real_t>
   int SYSWAPR(char uplo,int_t n,real_t *A,int_t ldA,int_t i,int_t j)
   {
      using std::toupper;

      uplo=toupper(uplo);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -4;
      else if((i<0)||(i>=n))
         return -5;
      else if((j<0)||(j>=n))
         return -6;
      
      if(uplo=='U') // upper triangular
      {
         real_t *ai=A+i*ldA;
         real_t *aj=A+j*ldA;
         real_t *ak;
         real_t t;
         SWAP(i,ai,1,aj,1);
         t=ai[i];
         ai[i]=aj[j];
         aj[j]=t;
         ak=ai;
         for(int_t k=1;k<j-i;k++)
         {
            ak+=ldA;
            t=ak[i];
            ak[i]=aj[i+k];
            aj[i+k]=t;
         }
         ak=aj;
         for(int_t k=j+1;k<n;k++)
         {
            ak+=ldA;
            t=ak[i];
            ak[i]=ak[j];
            ak[j]=t;
         }
      }
      else // lower triangular
      {
         real_t *ai=A+i*ldA;
         real_t *aj=A+j*ldA;
         real_t *ak;
         real_t t;
         SWAP(i,A+i,ldA,A+j,ldA);
         t=ai[i];
         ai[i]=aj[j];
         aj[j]=t;
         ak=ai;
         for(int_t k=1;k<j-i;k++)
         {
            ak+=ldA;
            t=ai[i+k];
            ai[i+k]=ak[j];
            ak[j]=t;
         }
         for(int_t k=j+1;k<n;k++)
         {
            t=ai[k];
            ai[k]=aj[k];
            aj[k]=t;
         }
      }
      return 0;
   }
}

#endif
