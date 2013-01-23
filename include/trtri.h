//
//  tritri.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _tritri_h
#define _tritri_h

/// @file trtri.h Computes the inverse of an upper or lower triangular matrix using blocked algorithm.

#include "latl.h"
#include "trti2.h"
#include "trmm.h"
#include "trsm.h"

namespace latl
{
   /// @brief Computes the inverse of an upper or lower triangular matrix using blocked algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i if the ith diagonal element is exactly zero,
   /// thus the matrix is singular and its inverse can not be computed.
   /// @param uplo Specifies whether the triangular factor stored in the array 
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param diag Specifies whether the or not the matrix is unit triangular:
   ///
   ///             'N' or 'n':  non-unit triangular
   ///             'U' or 'u':  unit triangular
   /// @param n Order of the triangular matrix A.  n >= 0.
   /// @param A Real triangular matrix of order n.
   /// If A is upper triangular, the strictly lower triangular part of A is not referenced.
   /// If A is lower triangular, the strictly upper triangular part of A is not referenced.
   /// If A is unit triangular, the diagonal elements of A are also not referenced and
   /// are assumed to be one. On exit, A contains the inverse of the original matrix.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param nb Block size (optional).
   /// @ingroup MATM

   template <typename real_t>
   int_t trtri(char uplo, char diag, int_t n, real_t *A, int_t ldA, int_t nb=32)
   {
      using std::toupper;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const real_t one=1.0;
      const real_t zero=0.0;
      int_t info;
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      if((diag!='N')&&(diag!='U'))
         return -2;
      else if(n<0)
         return -3;
      else if(ldA<n)
         return -5;
      else if(n==0)
         return 0;

      if(diag=='N')
      {
         for(int_t info=0; info<n; info++)
            if(A[info+info*ldA]==zero)
               return info;
         info = 0;
      }
      if(nb>n)
         info=trti2( uplo, diag, n, A, ldA);
      else
      {
         if(uplo=='U')
         {
            for(int_t j=0; j< n; j+=nb)
            {
               int_t jb = std::min(nb,n-j);
               trmm('L', 'U', 'N', diag, j, jb, one, A, ldA, A+j*ldA, ldA);
               trsm('R', 'U', 'N', diag, j, jb, -one, A+j+j*ldA, ldA, A+j*ldA, ldA);
               info=trti2('U', diag, jb, A+j+j*ldA, ldA);
            }
         }
         else
         {
            int_t nn=((n-1)/nb)*nb;
            for(int_t j=nn;j>=0;j-=nb)
            {
               int_t jb = std::min(nb,n-j);
               if(j+jb < n)
               {
                  trmm('L', 'L', 'N', diag, n-j-jb, jb, one, A+(j+jb)+(j+jb)*ldA, ldA, A+(j+jb)+j*ldA, ldA);
                  trsm('R', 'L', 'N', diag, n-j-jb, jb, -one, A+j+j*ldA, ldA, A+(j+jb)+j*ldA, ldA);
               }
               info=trti2('L', diag, jb, A+j+j*ldA, ldA);
            }
         }
      }
      return info;
   }
   
   
   /// @brief Computes the inverse of an upper or lower triangular matrix using blocked algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i if the ith diagonal element is exactly zero,
   /// thus the matrix is singular and its inverse can not be computed.
   /// @param uplo Specifies whether the triangular factor stored in the array
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param diag Specifies whether the or not the matrix is unit triangular:
   ///
   ///             'N' or 'n':  non-unit triangular
   ///             'U' or 'u':  unit triangular
   /// @param n Order of the triangular matrix A.  n >= 0.
   /// @param A Complex triangular matrix of order n.
   /// If A is upper triangular, the strictly lower triangular part of A is not referenced.
   /// If A is lower triangular, the strictly upper triangular part of A is not referenced.
   /// If A is unit triangular, the diagonal elements of A are also not referenced and
   /// are assumed to be one. On exit, A contains the inverse of the original matrix.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param nb Block size (optional).
   /// @ingroup MATM

   template <typename real_t>
   int_t trtri(char uplo, char diag, int_t n, complex<real_t> *A, int_t ldA, int_t nb=32)
   {
      using std::toupper;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const complex<real_t> one=1.0;
      const complex<real_t> zero=0.0;
      int_t info;
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      if((diag!='N')&&(diag!='U'))
         return -2;
      else if(n<0)
         return -3;
      else if(ldA<n)
         return -5;
      else if(n==0)
         return 0;
      
      if(diag=='N')
      {
         for(int_t info=0; info<n; info++)
            if(A[info+info*ldA]==zero)
               return info;
         info = 0;
      }
      if(nb>n)
         info=trti2( uplo, diag, n, A, ldA);
      else
      {
         if(uplo=='U')
         {
            for(int_t j=0; j< n; j+=nb)
            {
               int_t jb = std::min(nb,n-j);
               trmm('L', 'U', 'N', diag, j, jb, one, A, ldA, A+j*ldA, ldA);
               trsm('R', 'U', 'N', diag, j, jb, -one, A+j+j*ldA, ldA, A+j*ldA, ldA);
               info=trti2('U', diag, jb, A+j+j*ldA, ldA);
            }
         }
         else
         {
            int_t nn=((n-1)/nb)*nb;
            for(int_t j=nn;j>=0;j-=nb)
            {
               int_t jb = std::min(nb,n-j);
               if(j+jb < n)
               {
                  trmm('L', 'L', 'N', diag, n-j-jb, jb, one, A+(j+jb)+(j+jb)*ldA, ldA, A+(j+jb)+j*ldA, ldA);
                  trsm('R', 'L', 'N', diag, n-j-jb, jb, -one, A+j+j*ldA, ldA, A+(j+jb)+j*ldA, ldA);
               }
               info=trti2('L', diag, jb, A+j+j*ldA, ldA);
            }
         }
      }
      return info;
   }
}

#endif
