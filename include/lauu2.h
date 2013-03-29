//
//  lauu2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/17/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _lauu2_h
#define _lauu2_h

/// @file lauu2.h Computes product of a triangular matrix with its transpose using unblocked algorithm.

#include <cctype>
#include <algorithm>
#include "gemv.h"
#include "scal.h"
#include "dot.h"
#include "lacgv.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes product of a real triangular matrix with its transpose.
   ///
   /// The product U*U' or L'*L is computed, where the triangular factor U or L is stored in
   /// the upper or lower triangular part of the matrix A.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param uplo Specifies whether the triangular factor stored in the array is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param n The order of the triangular factor U or L.  n >= 0.
   /// @param A Real triangular matrix of order n.
   /// On entry, the triangular factor U or L.  On exit, if upper trianglar, A is overwritten with
   /// the upper triangle of the product U*U'; if lower trianglar, A is overwritten with the lower
   //  triangle of the product L'*L.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @ingroup AUX
   
   template <typename real_t>
   int LAUU2(char uplo, int_t n, real_t *A, int_t ldA)
   {
      using std::toupper;
      const real_t one(1.0);
      uplo=toupper(uplo);

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -4;
      else if(n==0)
         return 0;

      if(uplo=='U')
      {
         real_t *B=A+ldA;
         for(int_t i=0;i<n-1;i++)
         {
            real_t a=A[i];
            A[i]=DOT<real_t>(n-i,A+i,ldA,A+i,ldA);
            GEMV<real_t>('N',i,n-i-1,one,B,ldA,B+i,ldA,a,A,1);
            A+=ldA;
            B+=ldA;
         }
         SCAL<real_t>(n,A[n-1],A,1);
      }
      else
      {
         real_t *B=A;
         for(int_t i=0;i<n-1;i++)
         {
            real_t a=A[i];
            A[i]=DOT<real_t>(n-i,A+i,1,A+i,1);
            GEMV<real_t>('T',n-i-1,i,one,B+i+1,ldA,A+i+1,1,a,B+i,ldA);
            A+=ldA;
         }
         SCAL<real_t>(n,A[n-1],B+n-1,ldA);
      }
      return 0;
   }

   /// @brief Computes product of a complex triangular matrix with its conjugate transpose.
   ///
   /// The product U*U' or L'*L is computed, where the triangular factor U or L is stored in
   /// the upper or lower triangular part of the matrix A.
   ///
   ///        Note that the diagonal must be real.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param uplo Specifies whether the triangular factor stored in the array is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param n The order of the triangular factor U or L.  n >= 0.
   /// @param A Complex triangular matrix of order n, with real diagonal elements.
   /// On entry, the triangular factor U or L.  On exit, if upper trianglar, A is overwritten with
   /// the upper triangle of the product U*U'; if lower trianglar, A is overwritten with the lower
   //  triangle of the product L'*L.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @ingroup AUX

   template <typename real_t>
   int LAUU2(char uplo, int_t n, complex<real_t> *A, int_t ldA)
   {
      using std::toupper;
      using std::real;
      const complex<real_t> one(1.0);
      uplo=toupper(uplo);

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -4;
      else if(n==0)
         return 0;

      if(uplo=='U')
      {
         complex<real_t> *B=A+ldA;
         for(int_t i=0;i<n-1;i++)
         {
            real_t a=real(A[i]);
            A[i]=a*a+real(DOTC<real_t>(n-i-1,B+i,ldA,B+i,ldA));
            LACGV<real_t>(n-i-1,B+i,ldA);
            GEMV<real_t>('N',i,n-i-1,one,B,ldA,B+i,ldA,a,A,1);
            LACGV<real_t>(n-i-1,B+i,ldA);
            A+=ldA;
            B+=ldA;
         }
         real_t a=real(A[n-1]);
         SCAL<real_t>(n,a,A,1);
         A[n-1]=a*a;
      }
      else
      {
         complex<real_t> *B=A;
         for(int_t i=0;i<n-1;i++)
         {
            real_t a=real(A[i]);
            A[i]=a*a+real(DOTC<real_t>(n-i-1,A+i+1,1,A+i+1,1));
            LACGV<real_t>(i,B+i,ldA);
            GEMV<real_t>('C',n-i-1,i,one,B+i+1,ldA,A+i+1,1,a,B+i,ldA);
            LACGV<real_t>(i,B+i,ldA);
            A+=ldA;
         }
         real_t a=real(A[n-1]);
         SCAL<real_t>(n,a,B+n-1,ldA);
         A[n-1]=a*a;
      }
      return 0;
   }
}
#endif
