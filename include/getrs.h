//
//  getrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 5/30/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _getrs_h
#define _getrs_h

/// @file getrs.h Solves a system of linear equations A * x = B.

#include "laswp.h"
#include "trsm.h"
#include "latl.h"

namespace LATL
{
   
   /// @brief Solves a system of linear equations A * X = B.
   ///
   ///     A * X = B or A^T * X = B
   ///
   /// with a general n-by-n matrix A using the LU factorization computed by getrf.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies the form of the system of equations.  'N' for no transpose, 'T' for transpose, 'C' for Conjugate transpose.
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param colB Number of columns of the matrix B.  colB >= 0
   /// @param A Real matrix size ldA-by-n.  Should already contain factors L and U from the factorization A = P*L*U computed by getrf.
   /// @param ldA Column length of matrix A.  ldA >= n
   /// @param IPIV Permutation matrix size min(m,n), computed by getrf.  In A, row k of A was exchanged with IPIV[k].
   /// @param B Real matrix size ldB-by-colB.  On exit, the solution matrix X.
   /// @param ldB Column length of B.  ldB >= n
   /// @ingroup SOLV
   
   template<typename real_t>
   int GETRS( const char trans, const int_t n, const int_t colB, real_t * const A, const int_t ldA, int_t * const IPIV, real_t * const B, const int_t ldB)
   {
      bool notrans = ((trans == 'N') || (trans == 'n'));
      if ( !notrans && (trans != 'T') && (trans != 't') && (trans != 'C') && (trans != 'c'))
         return -1;
      if (n < 0)
         return -2;
      if (colB < 0)
         return -3;
      if (ldA < n)
         return -5;
      if (ldB < n)
         return -8;
      
      if ( n == 0 || colB == 0)
         return 0;
      
      const real_t one(1.0);
      
      if (notrans)
      {
         LATL::LASWP(colB, B, ldB, 0, n-1, IPIV);
         LATL::TRSM('L', 'L', 'N', 'U', n, colB, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'U', 'N', 'N', n, colB, one, A, ldA, B, ldB);
      }
      else
      {
         LATL::TRSM('L', 'U', 'T', 'N', n, colB, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'L', 'T', 'U', n, colB, one, A, ldA, B, ldB);
         LATL::LASWP(colB, B, ldB, 0, n-1, IPIV, -1);
      }
      return 0;
   }
   
   /// @brief Solves a system of linear equations A * X = B.
   ///
   ///     A * X = B or A^T * X = B
   ///
   /// with a general n-by-n matrix A using the LU factorization computed by getrf.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies the form of the system of equations.  'N' for no transpose, 'T' for transpose, 'C' for Conjugate transpose.
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param colB Number of columns of the matrix B.  colB >= 0
   /// @param A Complex matrix size ldA-by-n.  Should already contain factors L and U from the factorization A = P*L*U computed by getrf.
   /// @param ldA Column length of matrix A.  ldA >= n
   /// @param IPIV Permutation matrix size min(m,n), computed by getrf.  In A, row k of A was exchanged with IPIV[k].
   /// @param B Complex matrix size ldB-by-colB.  On exit, the solution matrix X.
   /// @param ldB Column length of B.  ldB >= n
   /// @ingroup SOLV

   template < typename real_t>
   int GETRS( const char trans, const int_t n, const int_t colB, complex<real_t> * const A, const int_t ldA, int_t * const IPIV, complex<real_t> * const B, const int_t ldB)
   {
      bool notrans = ((trans == 'N') || (trans == 'n'));
      if ( !notrans && (trans != 'T') && (trans != 't') && (trans != 'C') && (trans!= 'c'))
         return -1;
      if (n < 0)
         return -2;
      if (colB < 0)
         return -3;
      if (ldA < n)
         return -5;
      if (ldB < n)
         return -8;
      
      if ( n == 0 || colB == 0)
         return 0;
      
      const complex<real_t> one(1.0, 0);
      
      if (notrans)
      {
         LATL::LASWP(colB, B, ldB, 0, n-1, IPIV);
         LATL::TRSM('L', 'L', 'N', 'U', n, colB, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'U', 'N', 'N', n, colB, one, A, ldA, B, ldB);
      }
      else
      {
         LATL::TRSM('L', 'U', trans, 'N', n, colB, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'L', trans, 'U', n, colB, one, A, ldA, B, ldB);
         LATL::LASWP(colB, B, ldB, 0, n-1, IPIV, -1);
      }
      return 0;
   }
}

#endif
