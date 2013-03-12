//
//  gesv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 5/31/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _gesv_h
#define _gesv_h

/// @file gesv.h Computes the solution to a system of linear equations A * X = B.

#include "getrf.h"
#include "getrs.h"
#include "latl.h"

namespace latl
{
   /// @brief Computes the solution to a real system of linear equations
   ///
   ///     A * X = B
   ///
   /// where A is an n-by-n matrix and X and B are n-by-colB matrices.
   ///
   /// The LU decomposition with partial pivoting and row interchanges is used to factor A as
   ///
   ///     A = P * L * U
   ///
   /// where P is a permutation matrix, L is unit lower triangular, and U is upper triangular.  The factored form of A is then used to solve the system of equations.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column of matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param n Number of linear equations.  n >= 0
   /// @param colB Number of columns in matrix B.  colB >= 0
   /// @param A Real matrix size ldA-by-n.  On exit, the factors L and U.
   /// @param ldA Column length of the matrix A. ldA >= n
   /// @param IPIV Integer array of length n.  Defines the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).
   /// @param B Real matrix size n-by-colB.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B. ldB >= n
   /// @param nb Block size, optional.  Default value is 80.
   /// @ingroup SOLV
   
   template< typename real_t >
   int_t gesv(const int_t n, const int_t colB, real_t * const A, const int_t ldA, int_t * const IPIV, real_t * const B, const int_t ldB, int_t nb=32)
   {
      if (n < 0)
         return -1;
      if (colB < 0)
         return -2;
      if (ldA < n)
         return -4;
      if (ldB < n)
         return -7;
      
      int_t info = 0;
      info = latl::getrf<real_t>(n, n, A, ldA, IPIV, nb);
      if (info == 0)
      {
         info = latl::getrs('N', n, colB, A, ldA, IPIV, B, ldB);
      }
      return info;
   }
   
   /// @brief Computes the solution to a complex system of linear equations
   ///
   ///     A * X = B
   ///
   /// where A is an n-by-n matrix and X and B are n-by-colB matrices.
   ///
   /// The LU decomposition with partial pivoting and row interchanges is used to factor A as
   ///
   ///     A = P * L * U
   ///
   /// where P is a permutation matrix, L is unit lower triangular, and U is upper triangular.  The factored form of A is then used to solve the system of equations.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column of matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param n Number of linear equations.  n >= 0
   /// @param colB Number of columns in matrix B.  colB >= 0
   /// @param A Complex matrix size ldA-by-n.  On exit, the factors L and U.
   /// @param ldA Column length of the matrix A. ldA >= n
   /// @param IPIV Integer array of length n.  Defines the permutation matrix P; row i of the matrix was interchanged with row IPIV(i).
   /// @param B Complex matrix size n-by-colB.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B. ldB >= n
   /// @param nb Block size, optional.  Default value is 80.
   /// @ingroup SOLV
   
   template< typename real_t >
   int_t gesv(const int_t n, const int_t colB, complex<real_t> * const A, const int_t ldA, int_t * const IPIV, complex<real_t> * const B, const int_t ldB, int_t nb=32)
   {
      if (n < 0)
         return -1;
      if (colB < 0)
         return -2;
      if (ldA < n)
         return -4;
      if (ldB < n)
         return -7;
      
      int_t info = 0;
      info = latl::getrf<real_t>(n, n, A, ldA, IPIV, nb);
      if (info == 0)
      {
         info = latl::getrs('N', n, colB, A, ldA, IPIV, B, ldB);
      }
      return info;
   }
}

#endif
