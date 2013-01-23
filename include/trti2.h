//
//  triti2.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _triti2_h
#define _triti2_h

/// @file trti2.h Computes the inverse of a triangular matrix (unblocked algorithm).

#include "latl.h"
#include "trmv.h"
#include "scal.h"

namespace latl
{
   /// @brief Computes the inverse of a triangular matrix (unblocked algorithm).
   ///
   /// This is the Level 2 BLAS version of the algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -k if the k-th argument had an illegal value.
   /// @param uplo Specifies whether the triangular factor stored in the array 
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param diag Specifies whether the or not the matrix is unit triangular:
   ///
   ///             'N' or 'n':  non-unit triangular
   ///             'U' or 'u':  unit triangular
   /// @param n The order of the triangular factor U or L.  n >= 0.
   /// @param A Real triangular matrix of order n.
   /// A is DOUBLE PRECISION array, dimension (LDA,N)
   /// On entry, the triangular matrix A.  If UPLO = 'U', the leading n by n 
   /// upper triangular part of the array A contains the upper triangular 
   /// matrix, and the strictly lower triangular part of A is not referenced.
   /// If UPLO = 'L', the leading n by n lower triangular part of the array A 
   /// contains the lower triangular matrix, and the strictly upper triangular 
   /// part of A is not referenced.  If DIAG = 'U', the diagonal elements of A 
   /// are also not referenced and are assumed to be 1.
   ///
   /// On exit, the (triangular) inverse of the original matrix, in
   /// the same storage format.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @ingroup MATM

   template <typename real_t>
   int_t trti2(char uplo, char diag, int_t n, real_t *A, int_t ldA)
   {
      using std::toupper;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const real_t one=1.0;
      
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
      
      real_t AJJ;
      if(uplo=='U')
      {
         for(int_t j=0;j<n;j++)
         {
            if(diag=='N')
            {
               A[j+j*ldA] = one / A[j+j*ldA];
               AJJ = -A[j+j*ldA ];
            }
            else
               AJJ = -one;
            
            int_t info=trmv<real_t>('U','N', diag, j, A, ldA, A+j*ldA, 1);
            if(info!=0) return info;
            scal<real_t>(j,AJJ,A+j*ldA,1);
         }
      }
      else
      {
         for(int_t j=n-1;j>=0;j--)
         {
            if(diag=='N')
            {
               A[j+j*ldA] = one / A[j+j*ldA];
               AJJ = -A[j+j*ldA ];
            }
            else
               AJJ = -one;
            
            if(j<n-1)
            {
               int_t info=trmv('L', 'N', diag, n-j-1, A+(j+1)+(j+1)*ldA, ldA, A+(j+1)+j*ldA, 1);
               if(info!=0) return info;
               scal<real_t>(n-j-1, AJJ, A+(j+1)+j*ldA, 1);
            }
         }
      }
      return 0;
   }

   /// @brief Computes the inverse of a triangular matrix (unblocked algorithm).
   ///
   /// This is the Level 2 BLAS version of the algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -k if the k-th argument had an illegal value.
   /// @param uplo Specifies whether the triangular factor stored in the array 
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param diag Specifies whether the or not the matrix is unit triangular:
   ///
   ///             'N' or 'n':  non-unit triangular
   ///             'U' or 'u':  unit triangular
   /// @param n The order of the triangular factor U or L.  n >= 0.
   /// @param A Complex triangular matrix of order n.
   /// A is COMPLEX DOUBLE PRECISION array, dimension (LDA,N)
   /// On entry, the triangular matrix A.  If UPLO = 'U', the leading n by n 
   /// upper triangular part of the array A contains the upper triangular 
   /// matrix, and the strictly lower triangular part of A is not referenced.
   /// If UPLO = 'L', the leading n by n lower triangular part of the array A 
   /// contains the lower triangular matrix, and the strictly upper triangular 
   /// part of A is not referenced.  If DIAG = 'U', the diagonal elements of A 
   /// are also not referenced and are assumed to be 1.
   ///
   /// On exit, the (triangular) inverse of the original matrix, in
   /// the same storage format.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @ingroup MATM

   template <typename real_t>
   int_t trti2(char uplo, char diag, int_t n, complex<real_t> *A, int_t ldA)
   {
      using std::toupper;
      using std::real;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const complex<real_t> one=1.0;
      
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
      
      complex<real_t> AJJ;
      if(uplo=='U')
      {
         for(int_t j=0;j<n;j++)
         {
            if(diag=='N')
            {
               A[j+j*ldA] = one / A[j+j*ldA];
               AJJ = -A[j+j*ldA ];
            }
            else
               AJJ = -one;
            
            int_t info=trmv<real_t>('U','N', diag, j, A, ldA, A+j*ldA, 1);
            if(info!=0) return info;
            scal<real_t>(j,AJJ,A+j*ldA,1);
         }
      }
      else
      {
         for(int_t j=n-1;j>=0;j--)
         {
            if(diag=='N')
            {
               A[j+j*ldA] = one / A[j+j*ldA];
               AJJ = -A[j+j*ldA ];
            }
            else
               AJJ = -one;
            
            if(j<n-1)
            {
               int_t info=trmv('L', 'N', diag, n-j-1, A+(j+1)+(j+1)*ldA, ldA, A+(j+1)+j*ldA, 1);
               if(info!=0) return info;
               scal<real_t>(n-j-1, AJJ, A+(j+1)+j*ldA, 1);
            }
         }
      }
      return 0;
   }
}

#endif