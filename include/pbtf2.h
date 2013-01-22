//
//  pbtf2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/20/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _pbtf2_h
#define _pbtf2_h

/// @file pbtf2.h Computes the Cholesky factorization of a symmetric positive definite band matrix A.

#include "scal.h"
#include "syr.h"
#include "her.h"
#include "lacgv.h"
#include "latl.h"

namespace latl
{
   /// @brief Computes the Cholesky factorization of a real symmetric positive definite band matrix A.
   ///
   /// The factorization has the form
   ///
   ///         A = U' * U   if uplo = 'U'
   ///         A = L * L'   if uplo = 'L'
   ///
   /// where U is upper triangular, U' is the transpose of U, and L is lower triangular.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param kd The number of super-diagonals of the matrix A if uplo = 'U' or the number of sub-diagonals if uplo = 'L'.  kd >= 0
   /// @param AB Real array size ldAB-by-n.  On entry, the upper or lower triangle of the symmetric band matrix A, stored in the first kd+1 rows of the array. On exit, the triangular factor U or L from the Cholesky factorization A = U' * U or A = L * L' of the band matrix A in the same storage format as A.
   /// @param ldAB Column length of the array AB.  ldAB >= kd+1;
   /// @ingroup TRF
   
   template< typename real_t>
   int_t pbtf2(const char uplo, const int_t n, const int_t kd, real_t * const AB, const int_t ldAB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if ( n < 0)
         return -2;
      if (kd < 0 )
         return -3;
      if (ldAB < kd+1)
         return -5;
      
      if ( n == 0)
         return 0;
      
      const real_t one(1.0);
      const real_t zero(0.0);
      int_t kn, kld = ldAB-1;
      real_t ajj;
      real_t * ABj = AB;
      real_t * ABjp1 = AB+ldAB;
      
      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t j = 0; j < n; ++j)
         {
            ajj = ABj[kd];
            if (ajj <= zero || std::isnan(ajj))
               return j+1;
            ABj[kd] = ajj = sqrt(ajj);
            
            kn = std::min(kd, n-j-1);
            if (kn > 0)
            {
               latl::scal(kn, one/ajj, ABjp1+kd-1, kld);
               latl::syr('U', kn, -one, ABjp1+kd-1, kld, ABjp1+kd, kld);
            }
            ABj += ldAB;
            ABjp1 += ldAB;
         }
      }
      else
      {
         for (int_t j = 0; j < n; ++j)
         {
            ajj = ABj[0];
            if (ajj <= 0 || std::isnan(ajj))
               return j+1;
            ABj[0] = ajj = sqrt(ajj);
            kn = std::min(kd, n-j-1);
            if (kn > 0)
            {
               latl::scal(kn, one/ajj, ABj+1, 1);
               latl::syr('L', kn, -one, ABj+1, 1, ABjp1, kld);
            }
            ABj += ldAB;
            ABjp1 += ldAB;
         }
      }
      return 0;
   }
   
   /// @brief Computes the Cholesky factorization of a complex Hermitian positive definite band matrix A.
   ///
   /// The factorization has the form
   ///
   ///         A = U^H * U   if uplo = 'U'
   ///         A = L * L^H   if uplo = 'L'
   ///
   /// where U is upper triangular and U^H is the conjugate transpose of U, and L is lower triangular.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the Hermitian matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param kd The number of super-diagonals of the matrix A if uplo = 'U' or the number of sub-diagonals if uplo = 'L'.  kd >= 0
   /// @param AB Complex array size ldAB-by-n.  On entry, the upper or lower triangle of the Hermitian band matrix A, stored in the first kd+1 rows of the array.  On exit, the triangular factor U or L from the Cholesky factorization A = U^H * U or A = L * L^H of the band matrix A in the same storage format as A.
   /// @param ldAB Column length of the array AB.  ldAB >= kd+1;
   /// @ingroup TRF

   template< typename real_t>
   int_t pbtf2(const char uplo, const int_t n, const int_t kd, complex<real_t> * const AB, const int_t ldAB)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if ( n < 0)
         return -2;
      if (kd < 0 )
         return -3;
      if (ldAB < kd+1)
         return -5;
      
      if ( n == 0)
         return 0;
      
      const real_t one(1.0);
      const real_t zero(0.0);
      int_t kn, kld = ldAB-1;
      real_t ajj;
      complex<real_t> * ABj = AB;
      complex<real_t> * ABjp1 = AB+ldAB;
      
      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t j = 0; j < n; ++j)
         {
            ajj = real(ABj[kd]);
            if (ajj <= zero || std::isnan(ajj))
            {
               ABj[kd] = ajj;
               return j+1;
            }
            ABj[kd] = ajj = sqrt(ajj);
            
            kn = std::min(kd, n-j-1);
            if (kn > 0)
            {
               latl::scal(kn, one/ajj, ABjp1+kd-1, kld);
               latl::lacgv(kn, ABjp1+kd-1, kld);
               latl::her('U', kn, -one, ABjp1+kd-1, kld, ABjp1+kd, kld);
               latl::lacgv(kn, ABjp1+kd, kld);
            }
            ABj += ldAB;
            ABjp1 += ldAB;
         }
      }
      else
      {
         for (int_t j = 0; j < n; ++j)
         {
            ajj = real(ABj[0]);
            if (ajj <= 0 || std::isnan(ajj))
            {
               ABj[0] = ajj;
               return j+1;
            }
            ABj[0] = ajj = sqrt(ajj);
            kn = std::min(kd, n-j-1);
            if (kn > 0)
            {
               latl::scal(kn, one/ajj, ABj+1, 1);
               latl::her('L', kn, -one, ABj+1, 1, ABjp1, kld);
            }
            ABj += ldAB;
            ABjp1 += ldAB;
         }
      }
      return 0;
   }
}

#endif
