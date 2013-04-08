//
//  potrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/12/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _potrs_h
#define _potrs_h

/// @file potrs.h Solves a system of linear equations A * X = B.

#include "trsm.h"
#include "latl.h"

namespace LATL
{
   /// @brief Solves a real system of linear equations A * X = B
   ///
   ///   A * X = B
   ///
   /// where A is symmetric positive definite using the Cholesky factorization A = U' * U or A = L*L' computed by POTRF.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param A Real symmetric matrix size ldA-by-n.  On entry, the matrix A.  On exit, the factor U or L from the Cholesky factorization A = U' * U or A = L * L'.
   /// @param ldA Column length of matrix A. ldA >= n
   /// @param B Real matrix size n-by-colB.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B. ldB >= n
   /// @ingroup COMP
   
   template< typename real_t>
   int_t POTRS( const char uplo, const int_t n, const int_t nrhs, real_t * const A, const int_t ldA, real_t * const B, const int_t ldB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldA < n)
         return -5;
      if (ldB < n)
         return -7;
      
      if (n == 0 || nrhs == 0)
         return 0;
      const real_t one(1.0);
      if (uplo == 'U' || uplo == 'u')
      {
         LATL::TRSM('L', 'U', 'T', 'N', n, nrhs, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'U', 'N', 'N', n, nrhs, one, A, ldA, B, ldB);
      }
      else
      {
         LATL::TRSM('L', 'L', 'N', 'N', n, nrhs, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'L', 'T', 'N', n, nrhs, one, A, ldA, B, ldB);
      }
      return 0;
   }
   /// @brief Solves a complex system of linear equations A * X = B
   ///
   ///   A * X = B
   ///
   /// where A is symmetric positive definite using the Cholesky factorization A = U' * U or A = L*L' computed by POTRF.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param A Complex symmetric matrix size ldA-by-n.  On entry, the matrix A.  On exit, the factor U or L from the Cholesky factorization A = U' * U or A = L * L'.
   /// @param ldA Column length of matrix A. ldA >= n
   /// @param B Complex matrix size n-by-colB.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B. ldB >= n
   /// @ingroup COMP

   template< typename real_t>
   int_t POTRS( const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, const int_t ldA, complex<real_t> * const B, const int_t ldB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldA < n)
         return -5;
      if (ldB < n)
         return -7;
      
      if (n == 0 || nrhs == 0)
         return 0;
      const complex<real_t> one(1.0);
      if (uplo == 'U' || uplo == 'u')
      {
         LATL::TRSM('L', 'U', 'C', 'N', n, nrhs, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'U', 'N', 'N', n, nrhs, one, A, ldA, B, ldB);
      }
      else
      {
         LATL::TRSM('L', 'L', 'N', 'N', n, nrhs, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'L', 'C', 'N', n, nrhs, one, A, ldA, B, ldB);
      }
      return 0;
   }
}

#endif
