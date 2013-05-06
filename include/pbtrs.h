//
//  pbtrs.h
//  Linear Algebra Template Library  
//
//  Created by Stephanie Patterson on 4/29/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.//

#ifndef _pbtrs_h
#define _pbtrs_h

#include "tbsv.h"
#include "latl.h"

/// @file pbtrs.h  Computes the solution to a system of equations A * X = B.

namespace LATL
{
   /// @brief Computes the solution to a real system of equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n symmetric positive definite band matrix and X and B are n-by-nrhs matrices.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param kd The number of super-diagonals of the matrix A if uplo == 'U' or the number of sub-diagonals if uplo == 'L'.  kd >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param AB Real array size ldAB-by-n.  On entry, the Cholesky factorization of the symmetric band matrix A as computed by PBTRF.
   /// @param ldAB Column length of the matrix AB.  ldAB >= kd+1
   /// @param B Real array size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B. ldB >= n
   /// @ingroup COMP
   
   template< typename real_t>
   int_t PBTRS(const char uplo, const int_t n, const int_t kd, const int_t nrhs, real_t * const AB, const int_t ldAB, real_t * const B, const int_t ldB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if ( n < 0)
         return -2;
      if (kd < 0 )
         return -3;
      if (nrhs < 0)
         return -4;
      if (ldAB < kd+1)
         return -6;
      if (ldB < n)
         return -8;
      
      if ( n == 0 || nrhs == 0)
         return 0;
      
      real_t * Bj = B;
      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t j = 0; j < nrhs; ++j)
         {
            LATL::TBSV('U', 'T', 'N', n, kd, AB, ldAB, Bj, 1);
            LATL::TBSV('U', 'N', 'N', n, kd, AB, ldAB, Bj, 1);
            Bj += ldB;
         }
      }
      else
      {
         for (int_t j = 0; j < nrhs; ++j)
         {
            LATL::TBSV('L', 'N', 'N', n, kd, AB, ldAB, Bj, 1);
            LATL::TBSV('L', 'T', 'N', n, kd, AB, ldAB, Bj, 1);
            Bj += ldB;
         }
      }
      return 0;
   }
   
   /// @brief Computes the solution to a complex system of equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n symmetric positive definite band matrix and X and B are n-by-nrhs matrices.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param kd The number of super-diagonals of the matrix A if uplo == 'U' or the number of sub-diagonals if uplo == 'L'.  kd >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param AB Complex array size ldAB-by-n.  On entry, the Cholesky factorization of the symmetric band matrix A as computed by PBTRF.
   /// @param ldAB Column length of the matrix AB.  ldAB >= kd+1
   /// @param B Complex array size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B. ldB >= n
   /// @ingroup COMP
   
   template< typename real_t>
   int_t PBTRS(const char uplo, const int_t n, const int_t kd, const int_t nrhs, complex<real_t> * const AB, const int_t ldAB, complex<real_t> * const B, const int_t ldB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if ( n < 0)
         return -2;
      if (kd < 0 )
         return -3;
      if (nrhs < 0)
         return -4;
      if (ldAB < kd+1)
         return -6;
      if (ldB < n)
         return -8;
      
      if ( n == 0 || nrhs == 0)
         return 0;
      
      real_t * Bj = B;
      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t j = 0; j < nrhs; ++j)
         {
            LATL::TBSV('U', 'C', 'N', n, kd, AB, ldAB, Bj, 1);
            LATL::TBSV('U', 'N', 'N', n, kd, AB, ldAB, Bj, 1);
            Bj += ldB;
         }
      }
      else
      {
         for (int_t j = 0; j < nrhs; ++j)
         {
            LATL::TBSV('L', 'N', 'N', n, kd, AB, ldAB, Bj, 1);
            LATL::TBSV('L', 'C', 'N', n, kd, AB, ldAB, Bj, 1);
            Bj += ldB;
         }
      }
      return 0;
   }
   
}

#endif
