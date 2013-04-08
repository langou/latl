//
//  posv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/14/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _posv_h
#define _posv_h

/// @file posv.h  Computes the solution to a system of equations A * X = B.

#include "potrs.h"
#include "potrf.h"

namespace LATL
{
   /// @brief Computes the solution to a real system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n symmetric positive definite matrix and X and B are n-by-nrhs matrices.
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
   /// @param A Real symmetric matrix size ldA-by-n.  On entry, the matrix A.  On exit, the factor U or L from the Cholesky factorization A = U' * U or A = L * L'.
   /// @param ldA Column length of matrix A. ldA >= n
   /// @param B Real matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B.  ldB >= n
   /// @ingroup DRIV
   
   template<typename real_t>
   int_t POSV(const char uplo, const int_t n, const int_t nrhs, real_t * const A, const int_t ldA, real_t * const B, const int_t ldB)
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
      int_t info = 0;
      info = LATL::POTRF(uplo, n, A, ldA);
      if (info == 0)
          info = LATL::POTRS(uplo, n, nrhs, A, ldA, B, ldB);
      return info;
   }

   /// @brief Computes the solution to a complex system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n symmetric positive definite matrix and X and B are n-by-nrhs matrices.
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
   /// @param A Complex symmetric matrix size ldA-by-n.  On entry, the matrix A.  On exit, the factor U or L from the Cholesky factorization A = U' * U or A = L * L'.
   /// @param ldA Column length of matrix A. ldA >= n
   /// @param B Complex matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B.  ldB >= n
   /// @ingroup DRIV
   
   template<typename real_t>
   int_t POSV(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, const int_t ldA, complex<real_t> * const B, const int_t ldB)
   {
      return LATL::POSV< complex<real_t> >(uplo, n, nrhs, A, ldA, B, ldB);
   }
}

#endif
