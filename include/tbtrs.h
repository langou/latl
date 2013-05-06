//
//  tbtrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 5/6/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _tbtrs_h
#define _tbtrs_h

/// @file tbtrs.h  Computes the solution to a system of equations A * X = B.

#include "latl.h"
#include "tbsv.h"

namespace LATL
{
   /// @brief Computes the solution to a real system of equations
   ///
   ///      A * X = B      or       A' * X = B
   ///
   /// where A is an n-by-n triangular band matrix and X and B are n-by-nrhs matrices.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith pivot is 0.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the band matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param trans Specifies the form of the system of equations:
   ///         A * X = B      if trans == 'N'
   ///         A' * X = B     if trans == 'T' or 'C'
   /// @param diag Specifies if A is unit triangular.  
   /// @param n Order of the matrix A.  n >= 0
   /// @param kd The number of super-diagonals of the matrix A if uplo == 'U' or the number of sub-diagonals if uplo == 'L'.  kd >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param AB Real array size ldAB-by-n.  On entry, the triangular band matrix A stored in the first kd+1 rows of AB. 
   /// @param ldAB Column length of the matrix AB.  ldAB >= kd+1
   /// @param B Real array size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B. ldB >= n
   /// @ingroup COMP
   
   template< typename real_t>
   int_t TBTRS(const char uplo, const char trans, const char diag, const int_t n, const int_t kd, const int_t nrhs, real_t * const AB, const int_t ldAB, real_t * const B, const int_t ldB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (trans != 'N' && trans != 'T' && trans != 'C' && trans != 'n' && trans != 't' && trans != 'c')
         return -2;
      if ( diag != 'N' && diag != 'U' && diag != 'n' && diag != 'u')
         return -3;
      if (n < 0)
         return -4;
      if (kd < 0)
         return -5;
      if (nrhs < 0)
         return -6;
      if ( ldAB < kd+1)
         return -8;
      if (ldB < n)
         return -10;
      
      if ( n == 0)
         return 0;
      
      real_t * ABj = AB;
      const real_t zero(0.0);
      
      if (diag == 'N' || diag == 'n')
      {
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (ABj[kd] == zero)
                  return j+1;
               ABj += ldAB;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (ABj[0] == zero)
                  return j+1;
               ABj += ldAB;
            }
         }
      }
      
      real_t * Bj = B;
      for (int_t j = 0; j < nrhs; ++j)
      {
         LATL::TBSV(uplo, trans, diag, n, kd, AB, ldAB, Bj, 1);
         Bj += ldB;
      }
      return 0;
   }
   
   /// @brief Computes the solution to a complex system of equations
   ///
   ///      A * X = B      or       A' * X = B     or    A.' * X = B
   ///
   /// where A is an n-by-n triangular band matrix and X and B are n-by-nrhs matrices.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith pivot is 0.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the band matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param trans Specifies the form of the system of equations:
   ///         A * X = B      if trans == 'N'
   ///         A' * X = B     if trans == 'C'
   ///         A.' * X = B    if trans == 'T'
   /// @param diag Specifies if A is unit triangular.
   /// @param n Order of the matrix A.  n >= 0
   /// @param kd The number of super-diagonals of the matrix A if uplo == 'U' or the number of sub-diagonals if uplo == 'L'.  kd >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param AB Complex array size ldAB-by-n.  On entry, the triangular band matrix A stored in the first kd+1 rows of AB.
   /// @param ldAB Column length of the matrix AB.  ldAB >= kd+1
   /// @param B Complex array size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B. ldB >= n
   /// @ingroup COMP
   template< typename real_t>
   int_t TBTRS(const char uplo, const char trans, const char diag, const int_t n, const int_t kd, const int_t nrhs, complex<real_t> * const AB, const int_t ldAB, complex<real_t> * const B, const int_t ldB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (trans != 'N' && trans != 'T' && trans != 'C' && trans != 'n' && trans != 't' && trans != 'c')
         return -2;
      if ( diag != 'N' && diag != 'U' && diag != 'n' && diag != 'u')
         return -3;
      if (n < 0)
         return -4;
      if (kd < 0)
         return -5;
      if (nrhs < 0)
         return -6;
      if ( ldAB < kd+1)
         return -8;
      if (ldB < n)
         return -10;
      
      if ( n == 0)
         return 0;
      
      complex<real_t> * ABj = AB;
      const complex<real_t> zero(0.0);
      
      if (diag == 'N' || diag == 'n')
      {
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (ABj[kd] == zero)
                  return j+1;
               ABj += ldAB;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (ABj[0] == zero)
                  return j+1;
               ABj += ldAB;
            }
         }
      }
      
      complex<real_t> * Bj = B;
      for (int_t j = 0; j < nrhs; ++j)
      {
         LATL::TBSV(uplo, trans, diag, n, kd, AB, ldAB, Bj, 1);
         Bj += ldB;
      }
      return 0;
   }
}

#endif
