//
//  potri.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 1/23/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _potri_h
#define _potri_h

/// @file potri.h Computes the inverse of a symmetric positive definite matrix.

#include "latl.h"
#include "trtri.h"
#include "lauum.h"

namespace LATL
{
   /// @brief Computes the inverse of a symmetric positive definite matrix.
   /// 
   /// Using the Cholesky factorization of A as computed by POTRF, the inverse
   /// of the symmetric positive definite matrix is returned in either the 
   /// upper, U, or lower, L, part of the matrix A.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @param uplo Specifies whether the triangular factor stored in the array
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param n The order of the triangular factor U or L.  n >= 0.
   /// @param A Real triangular matrix of order n.
   /// On entry, the triangular Cholesky factorization U or L.  On exit, if 
   /// upper trianglar, A is overwritten with the upper triangle of the inverse
   /// of A; if lower trianglar, A is overwritten with the lower triangle of 
   /// the inverse of A.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param nb Block size (optional).
   /// @ingroup MATM

   template <typename real_t>
      int_t POTRI(char uplo, int_t n, real_t *A, int_t ldA, int_t nb=40)
      {

         using std::toupper;
         uplo=toupper(uplo);
         int_t info;

         if((uplo!='U')&&(uplo!='L'))
            return -1;
         else if(n<0)
            return -2;
         else if(ldA<n)
            return -4;
         else if(n==0)
            return 0;

         info=TRTRI<real_t>( uplo, 'N', n, A, ldA, nb);
         if (info != 0) 
            return info;
         info=LAUUM<real_t>( uplo, n, A, ldA, nb);
         return info;
      }

   /// @brief Computes the inverse of a Hermitian positive definite matrix.
   /// 
   /// Using the Cholesky factorization of A as computed by POTRF, the inverse
   /// of the Hermitian positive definite matrix is returned in either the 
   /// upper, U, or lower, L, part of the matrix A.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @param uplo Specifies whether the triangular factor stored in the array
   /// is upper or lower triangular:
   ///
   ///             'U' or 'u':  upper triangular
   ///             'L' or 'l':  lower triangular
   /// @param n The order of the triangular factor U or L.  n >= 0.
   /// @param A complex triangular matrix of order n.
   /// On entry, the triangular Cholesky factorization U or L.  On exit, if 
   /// upper trianglar, A is overwritten with the upper triangle of the inverse
   /// of A; if lower trianglar, A is overwritten with the lower triangle of 
   /// the inverse of A.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param nb Block size (optional).
   /// @ingroup MATM

   template <typename real_t>
      int_t POTRI(char uplo, int_t n, complex<real_t> *A, int_t ldA, int_t nb=40)
      {

         using std::toupper;
         uplo=toupper(uplo);
         int_t info;

         if((uplo!='U')&&(uplo!='L'))
            return -1;
         else if(n<0)
            return -2;
         else if(ldA<n)
            return -4;
         else if(n==0)
            return 0;

         info=TRTRI<real_t>( uplo, 'N', n, A, ldA, nb);
         if (info != 0) 
            return info;
         info=LAUUM<real_t>( uplo, n, A, ldA, nb);
         return info;
      }
}

#endif
