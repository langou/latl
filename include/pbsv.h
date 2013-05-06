//
//  pbsv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 5/6/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _pbsv_h
#define _pbsv_h

/// @file pbsv.h  Solves a system of linear equations A * X = B.

#include "latl.h"
#include "pbtrf.h"
#include "pbtrs.h"

namespace LATL
{
   /// @brief Computes the solution to a real system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n symmetric positive definite band matrix and X and B are n-by-nrhs matrices.
   ///
   /// The Cholesky decomposition is used to factor A as
   ///      A = U'*U    if uplo == 'U'
   ///      A = L * L'  if uplo == 'L'
   /// where U is an upper triangular matrix and L is a lower triangular matrix.  The factored form of A is then used to solve the system of equations A * X = B.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the matrix A is stored as upper triangular or lower triangular.
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of B.  nrhs >= 0
   /// @param AB Real array size ldAB-by-n.  On entry, the upper or lower triangle of the symmetric band matrix A sotred in the first kd+1 rows of the array.  On exit, the factor U or L from the Cholesky factorization A = U' * U or A = L * L' in the same storage format as A.
   /// @param ldAB Column length of array AB. ldA >= n
   /// @param B Real matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B.  ldB >= n
   /// @ingroup DRIV
   
   template<typename real_t>
   int_t PBSV(const char uplo, const int_t n, const int_t kd, const int_t nrhs, real_t * const AB, const int_t ldAB, real_t * const B, const int_t ldB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (kd < 0)
         return -3;
      if (nrhs < 0)
         return -4;
      if (ldAB < kd+1)
         return -6;
      if (ldB < n)
         return -8;
      
      int_t info = LATL::PBTRF(uplo, n, kd, AB, ldAB);
      if ( info == 0)
      {
         info = LATL::PBTRS(uplo, n, kd, nrhs, AB, ldAB, B, ldB);
      }
      
      return info;
   }
   
   /// @brief Computes the solution to a complex system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n symmetric positive definite band matrix and X and B are n-by-nrhs matrices.
   ///
   /// The Cholesky decomposition is used to factor A as
   ///      A = U'*U    if uplo == 'U'
   ///      A = L * L'  if uplo == 'L'
   /// where U is an upper triangular matrix and L is a lower triangular matrix.  The factored form of A is then used to solve the system of equations A * X = B.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the matrix A is stored as upper triangular or lower triangular.
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of B.  nrhs >= 0
   /// @param AB Complex array size ldAB-by-n.  On entry, the upper or lower triangle of the symmetric band matrix A sotred in the first kd+1 rows of the array.  On exit, the factor U or L from the Cholesky factorization A = U' * U or A = L * L' in the same storage format as A.
   /// @param ldAB Column length of array AB. ldA >= n
   /// @param B Complex matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B.  ldB >= n
   /// @ingroup DRIV
   
   template<typename real_t>
   int_t PBSV(const char uplo, const int_t n, const int_t kd, const int_t nrhs, complex<real_t> * const AB, const int_t ldAB, complex<real_t> * const B, const int_t ldB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (kd < 0)
         return -3;
      if (nrhs < 0)
         return -4;
      if (ldAB < kd+1)
         return -6;
      if (ldB < n)
         return -8;
      
      int_t info = LATL::PBTRF(uplo, n, kd, AB, ldAB);
      if ( info == 0)
      {
         info = LATL::PBTRS(uplo, n, kd, nrhs, AB, ldAB, B, ldB);
      }
      
      return info;
   }
}

#endif
