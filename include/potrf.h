//
//  potrf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/20/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _potrf_h
#define _potrf_h

/// @file potrf.h Computes the Cholesky factorization of a symmetric positive definite matrix A.

#include "potf2.h"
#include "syrk.h"
#include "gemm.h"
#include "trsm.h"
#include "herk.h"
#include "latl.h"

namespace latl
{
   /// @brief Computes the Cholesky factorization of a real symmetric positive definite matrix A.
   ///
   /// The factorization has the form
   ///
   ///         A = U' * U  if uplo = 'U'
   ///         A = L * L'  if uplo = 'L'
   ///
   /// where U is an upper triangular matrix and L is lower triangular.  This is the blocked version of the algorithm.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Real symmetric matrix size ldA-by-n.  On entry, the matrix A.  On exit, the factor U or L from the Cholesky factorization A = U' * U or A = L * L'.
   /// @param ldA Column length of matrix A. ldA >= n
   /// @param nb Block size, optional.  Default value is 64.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t potrf(const char uplo, const int_t n, real_t * const A, const int_t ldA, const int_t nb = 64)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      if (nb <= 1 || nb >= n)
         return latl::potf2(uplo, n, A, ldA);
      else
      {
         int_t info = 0, jb;
         const real_t one(1.0);
         real_t * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; j+=nb)
            {
               jb = std::min( nb, n-j);
               latl::syrk('U', 'T', jb, j, -one, Aj, ldA, one, Aj+j, ldA);
               info = latl::potf2('U', jb, Aj+j, ldA);
               if (info != 0)
                  return info+j;
               if (j+jb < n)
               {
                  latl::gemm('T', 'N', jb, n-j-jb, j, -one, Aj, ldA, Aj+ldA*jb, ldA, one, Aj+ldA*jb+j, ldA);
                  latl::trsm('L', 'U', 'T', 'N', jb, n-j-jb, one, Aj+j, ldA, Aj+j+ldA*jb, ldA);
               }
               Aj += ldA*nb;
            }
         }
         else
         {
            for (int_t j = 0; j < n; j += nb)
            {
               jb = std::min(nb, n-j);
               latl::syrk('L', 'N', jb, j, -one, A+j, ldA, one, Aj+j, ldA);
               info = latl::potf2('L', jb, Aj+j, ldA);
               if (info != 0)
                  return info+j;
               if (j+jb < n)
               {
                  latl::gemm('N', 'T', n-j-jb, jb, j, -one, A+j+jb, ldA, A+j, ldA, one, Aj+j+jb, ldA );
                  latl::trsm('R', 'L', 'T', 'N', n-j-jb, jb, one, Aj+j, ldA, Aj+j+jb, ldA);
               }
               Aj += ldA*nb;
            }
         }
         return 0;
      }
   }
   
   /// @brief Computes the Cholesky factorization of a complex Hermitian positive definite matrix A.
   ///
   /// The factorization has the form
   ///
   ///     A = U^H * U if uplo = 'U'
   ///     A = L * L^H if uplo = 'L'
   ///
   /// where U is an upper triangular matrix and L is lower triangular.  This is the blocked version of the algorithm.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Hermitian matrix size ldA-by-n.  On entry, the matrix A.  On exit, the factor U or L from the Cholesky factorization A = U^H * U or A = L * L^H.
   /// @param ldA Column length of matrix A. ldA >= n
   /// @param nb Block size, optional.  Default value is 64.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t potrf(const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA, const int_t nb = 64)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      if (nb <= 1 || nb >= n)
         return latl::potf2(uplo, n, A, ldA);
      else
      {
         int_t info = 0, jb;
         const real_t one(1.0);
         const complex<real_t> onec(1.0);
         complex<real_t> * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; j+=nb)
            {
               jb = std::min( nb, n-j);
               latl::herk('U', 'C', jb, j, -one, Aj, ldA, one, Aj+j, ldA);
               info = latl::potf2('U', jb, Aj+j, ldA);
               if (info != 0)
                  return info+j;
               if (j+jb < n)
               {
                  latl::gemm('C', 'N', jb, n-j-jb, j, -onec, Aj, ldA, Aj+ldA*jb, ldA, onec, Aj+ldA*jb+j, ldA);
                  latl::trsm('L', 'U', 'C', 'N', jb, n-j-jb, onec, Aj+j, ldA, Aj+j+ldA*jb, ldA);
               }
               Aj += ldA*nb;
            }
         }
         else
         {
            for (int_t j = 0; j < n; j += nb)
            {
               jb = std::min(nb, n-j);
               latl::herk('L', 'N', jb, j, -one, A+j, ldA, one, Aj+j, ldA);
               info = latl::potf2('L', jb, Aj+j, ldA);
               if (info != 0)
                  return info+j;
               if (j+jb < n)
               {
                  latl::gemm('N', 'C', n-j-jb, jb, j, -onec, A+j+jb, ldA, A+j, ldA, onec, Aj+j+jb, ldA );
                  latl::trsm('R', 'L', 'C', 'N', n-j-jb, jb, onec, Aj+j, ldA, Aj+j+jb, ldA);
               }
               Aj += ldA*nb;
            }
         }
         return 0;
      }
   }
}

#endif
