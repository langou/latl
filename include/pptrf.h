//
//  pptrf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/22/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _pptrf_h
#define _pptrf_h

/// @file pptrf.h Computes the Cholesky factorization of a symmetric positive definite matrix A stored in packed form.

#include "hpr.h"
#include "scal.h"
#include "spr.h"
#include "tpsv.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes the Cholesky factorization of a real symmetric positive definite matrix A stored in packed form.
   ///
   /// The factorization has the form
   ///
   ///         A = L * L'  if uplo = 'L'
   ///         A = U' * U  if uplo = 'U'
   ///
   /// where U is an upper triangular matrix and L is lower triangular.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param AP Real array, size n-by-(n+1)/2.  On entry, the upper or lower triangle of the symmetric matrix A, packed columnwise in a linear array.  The jth column of A is stored in the array AP as:
   ///
   ///         AP[i + (j+1)*j/2] = A[i, j]      for 0 <= i <= j    if uplo = 'U'
   ///         AP[i + (j)*(2n-j-1)/2] = A[i, j] for j <= i <= n-1  if uplo = 'L'
   ///
   /// On exit, if return value equals 0, the triangular factor U or L from the Cholesky factorization in the same storage format as A.
   /// @ingroup COMP
   
   template< typename real_t>
   int_t PPTRF( const char uplo, const int_t n, real_t * const AP)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if ( n < 0)
         return -2;
      
      if ( n == 0)
         return 0;
      
      const real_t zero(0.0);
      int_t jj = 0;
      real_t ajj;
      
      if (uplo == 'U' || uplo == 'u')
      {
         
         int_t jc;
         jj = -1;
         for (int_t j = 0; j < n; ++j)
         {
            jc = jj+1;
            jj += (j+1);
            
            if ( j > 0)
               LATL::TPSV('U', 'T', 'N', j, AP, AP+jc, 1);
            ajj = AP[jj] - LATL::DOT(j, AP+jc, 1, AP+jc, 1);
            if (ajj <= zero)
            {
               AP[jj] = ajj;
               return j+1;
            }
            AP[jj] = sqrt(ajj);
         }
      }
      else
      {
         const real_t one(1.0);
         
         for (int_t j = 0; j < n; ++j)
         {
            ajj = AP[jj];
            if (ajj <= zero)
            {
               return j+1;
            }
            AP[jj] = ajj = std::sqrt(ajj);
            
            if (j < n)
            {
               LATL::SCAL(n-j-1, one/ajj, AP+jj+1, 1);
               LATL::SPR('L', n-j-1, -one, AP+jj+1, 1, AP+jj+n-j);
               jj = jj+n-j;
            }
         }
      }
      
      return 0;
   }
   
   /// @brief Computes the Cholesky factorization of a complex Hermitian positive definite matrix A stored in packed form.
   ///
   /// The factorization has the form
   ///
   ///         A = L * L'  if uplo = 'L'
   ///         A = U' * U  if uplo = 'U'
   ///
   /// where U is an upper triangular matrix and L is lower triangular.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param AP Real array, size n-by-(n+1)/2.  On entry, the upper or lower triangle of the Hermitian matrix A, packed columnwise in a linear array.  The jth column of A is stored in the array AP as:
   ///
   ///         AP[i + (j+1)*j/2] = A[i, j]      for 0 <= i <= j    if uplo = 'U'
   ///         AP[i + (j)*(2n-j-1)/2] = A[i, j] for j <= i <= n-1  if uplo = 'L'
   ///
   /// On exit, if return value equals 0, the triangular factor U or L from the Cholesky factorization in the same storage format as A.
   /// @ingroup COMP
   
   template< typename real_t>
   int_t PPTRF( const char uplo, const int_t n, complex<real_t> * const AP)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if ( n < 0)
         return -2;
      
      if ( n == 0)
         return 0;
      
      const real_t zero(0.0);
      int_t jj = 0;
      real_t ajj;
      
      if (uplo == 'U' || uplo == 'u')
      {
         
         int_t jc;
         jj = -1;
         for (int_t j = 0; j < n; ++j)
         {
            jc = jj+1;
            jj += (j+1);
            
            if ( j > 0)
               LATL::TPSV('U', 'C', 'N', j, AP, AP+jc, 1);
            ajj = real(AP[jj] - LATL::DOTC(j, AP+jc, 1, AP+jc, 1));
            if (ajj <= zero)
            {
               AP[jj] = ajj;
               return j+1;
            }
            AP[jj] = sqrt(ajj);
         }
      }
      else
      {
         const real_t one(1.0);
         
         for (int_t j = 0; j < n; ++j)
         {
            ajj = real(AP[jj]);
            if (ajj <= zero)
            {
               return j+1;
            }
            AP[jj] = ajj = std::sqrt(ajj);
            
            if (j < n)
            {
               LATL::SCAL(n-j-1, one/ajj, AP+jj+1, 1);
               LATL::HPR('L', n-j-1, -one, AP+jj+1, 1, AP+jj+n-j);
               jj = jj+n-j;
            }
         }
      }
      
      return 0;
   }
}

#endif
