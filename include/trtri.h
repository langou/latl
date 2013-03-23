//
//  trtri.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 1/22/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _trtri_h
#define _trtri_h

/// @file trtri.h Computes the inverse of an upper or lower triangular matrix.

#include "latl.h"
#include "trmm.h"
#include "trsm.h"
#include "trmv.h"
#include "scal.h"

namespace LATL
{
   /// @brief Computes the inverse of a triangular matrix.
   ///
   /// This is the unblocked version of the algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument had an illegal value.
   /// @param uplo Specifies whether the triangular factor stored in the array
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param diag Specifies whether or not the matrix is unit triangular:
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

   template <typename real_t>
   int_t TRTRI(char uplo, char diag, int_t n, real_t *A, int_t ldA)
   {
      using std::toupper;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const real_t one(1.0);

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

            TRMV<real_t>('U','N', diag, j, A, ldA, A+j*ldA, 1);
            SCAL<real_t>(j,AJJ,A+j*ldA,1);
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
               TRMV('L', 'N', diag, n-j-1, A+(j+1)+(j+1)*ldA, ldA, A+(j+1)+j*ldA, 1);
               SCAL<real_t>(n-j-1, AJJ, A+(j+1)+j*ldA, 1);
            }
         }
      }
      return 0;
   }

   /// @brief Computes the inverse of a triangular matrix.
   ///
   /// This is the unblocked version of the algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument had an illegal value.
   /// @param uplo Specifies whether the triangular factor stored in the array
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param diag Specifies whether or not the matrix is unit triangular:
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

   template <typename real_t>
   int_t TRTRI(char uplo, char diag, int_t n, complex<real_t> *A, int_t ldA)
   {
      using std::toupper;
      using std::real;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const complex<real_t> one(1.0);

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

            TRMV<real_t>('U','N', diag, j, A, ldA, A+j*ldA, 1);
            SCAL<real_t>(j,AJJ,A+j*ldA,1);
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
               TRMV('L', 'N', diag, n-j-1, A+(j+1)+(j+1)*ldA, ldA, A+(j+1)+j*ldA, 1);
               SCAL<real_t>(n-j-1, AJJ, A+(j+1)+j*ldA, 1);
            }
         }
      }
      return 0;
   }

   /// @brief Computes the inverse of an upper or lower triangular matrix using blocked algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return 1 if the matrix is singular.
   /// @param uplo Specifies whether the triangular factor stored in the array 
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param diag Specifies whether or not the matrix is unit triangular:
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
   /// @param nb Block size.

   template <typename real_t>
   int_t TRTRI(char uplo, char diag, int_t n, real_t *A, int_t ldA, int_t nb)
   {
      using std::toupper;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const real_t one(1.0);
      const real_t zero(0.0);
      
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
         for(int_t i=0; i<n; i++)
            if(A[i+i*ldA]==zero)
               return 1;
      }
      if((nb<2)||(nb>n))
         TRTRI<real_t>( uplo, diag, n, A, ldA);
      else
      {
         if(uplo=='U')
         {
            for(int_t j=0; j< n; j+=nb)
            {
               int_t jb = std::min(nb,n-j);
               TRMM<real_t>('L', 'U', 'N', diag, j, jb, one, A, ldA, A+j*ldA, ldA);
               TRSM<real_t>('R', 'U', 'N', diag, j, jb, -one, A+j+j*ldA, ldA, A+j*ldA, ldA);
               TRTRI<real_t>('U', diag, jb, A+j+j*ldA, ldA);
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
                  TRMM<real_t>('L', 'L', 'N', diag, n-j-jb, jb, one, A+(j+jb)+(j+jb)*ldA, ldA, A+(j+jb)+j*ldA, ldA);
                  TRSM<real_t>('R', 'L', 'N', diag, n-j-jb, jb, -one, A+j+j*ldA, ldA, A+(j+jb)+j*ldA, ldA);
               }
               TRTRI<real_t>('L', diag, jb, A+j+j*ldA, ldA);
            }
         }
      }
      return 0;
   }
   
   /// @brief Computes the inverse of an upper or lower triangular matrix using blocked algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return 1 if the matrix is singular.
   /// @param uplo Specifies whether the triangular factor stored in the array
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param diag Specifies whether or not the matrix is unit triangular:
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
   /// @param nb Block size.

   template <typename real_t>
   int_t TRTRI(char uplo, char diag, int_t n, complex<real_t> *A, int_t ldA, int_t nb)
   {
      using std::toupper;
      uplo=toupper(uplo);
      diag=toupper(diag);
      const complex<real_t> one(1.0);
      const complex<real_t> zero(0.0);
      
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
         for(int_t i=0; i<n; i++)
            if(A[i+i*ldA]==zero)
               return 1;
      }
      if((nb>=n)||(nb<2))
      {
         TRTRI<real_t>( uplo, diag, n, A, ldA);
      }
      else
      {
         if(uplo=='U')
         {
            for(int_t j=0; j< n; j+=nb)
            {
               int_t jb = std::min(nb,n-j);
               TRMM<real_t>('L', 'U', 'N', diag, j, jb, one, A, ldA, A+j*ldA, ldA);
               TRSM<real_t>('R', 'U', 'N', diag, j, jb, -one, A+j+j*ldA, ldA, A+j*ldA, ldA);
               TRTRI<real_t>('U', diag, jb, A+j+j*ldA, ldA);
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
                  TRMM<real_t>('L', 'L', 'N', diag, n-j-jb, jb, one, A+(j+jb)+(j+jb)*ldA, ldA, A+(j+jb)+j*ldA, ldA);
                  TRSM<real_t>('R', 'L', 'N', diag, n-j-jb, jb, -one, A+j+j*ldA, ldA, A+(j+jb)+j*ldA, ldA);
               }
               TRTRI<real_t>('L', diag, jb, A+j+j*ldA, ldA);
            }
         }
      }
      return 0;
   }
}

#endif
