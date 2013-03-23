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

#include "dot.h"
#include "gemv.h"
#include "scal.h"
#include "lacgv.h"
#include "syrk.h"
#include "gemm.h"
#include "trsm.h"
#include "herk.h"
#include "latl.h"
#include <cmath>

namespace LATL
{
   /// @brief Computes the Cholesky factorization of a real symmetric positive definite matrix A.
   ///
   /// The factorization has the form
   ///
   ///         A = U' * U  if uplo = 'U'
   ///         A = L * L'  if uplo = 'L'
   ///
   /// where U is an upper triangular matrix and L is lower triangular.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Real symmetric matrix size ldA-by-n.  On entry, the matrix A.  On exit, the factor U or L from the Cholesky factorization A = U' * U or A = L * L'.
   /// @param ldA Column length of matrix A. ldA >= n

   template< typename real_t>
   int_t potrf(const char uplo, const int_t n, real_t * A, const int_t ldA)
   {
      using std::isnan;
      using std::sqrt;
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;

      if (n == 0)
         return 0;

      const real_t one(1.0);
      const real_t zero(0.0);
      real_t ajj;
      real_t * Aj = A;
      int_t info = 0;
      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t j = 0; j < n; ++j)
         {
            ajj = Aj[j] - LATL::dot(j, Aj, 1, Aj, 1);
            if (ajj <= zero or isnan(ajj))
            {
               Aj[j] = ajj;
               return j+1;
            }
            Aj[j] = ajj = sqrt(ajj);

            if (j < n-1)
            {
               LATL::gemv('T', j, n-j-1, -one, Aj+ldA, ldA, Aj, 1, one, Aj+ldA+j, ldA);
               LATL::scal(n-j-1, one/ajj, Aj+ldA+j, ldA);
            }
            Aj += ldA;
         }
      }
      else
      {
         for (int_t j = 0; j < n; ++j)
         {
            ajj = Aj[j] - LATL::dot(j, A+j, ldA, A+j, ldA);
            if (ajj <= zero || isnan(ajj))
            {
               Aj[j] = ajj;
               return j+1;
            }
            Aj[j] = ajj = sqrt(ajj);

            if (j < n-1)
            {
               LATL::gemv('N', n-j-1, j, -one, A+j+1, ldA, A+j, ldA, one, Aj+j+1, 1);
               LATL::scal(n-j-1, one/ajj, Aj+j+1, 1);
            }
            Aj += ldA;
         }
      }
      return info;
   }

   /// @brief Computes the Cholesky factorization of a complex Hermitian positive definite matrix A.
   ///
   /// The factorization has the form
   ///
   ///         A = U^H * U if uplo = 'U'
   ///         A = L * L^H if uplo = 'L'
   ///
   /// where U is an upper triangular matrix and L is lower triangular.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the Hermitian matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Hermitian matrix size ldA-by-n.  On entry, the matrix A.  On exit, the factor U or L from the Cholesky factorization A = U^H * U or A = L * L^H.
   /// @param ldA Column length of matrix A. ldA >= n

   template< typename real_t>
   int_t potrf(const char uplo, const int_t n, complex<real_t> * A, const int_t ldA)
   {
      using std::sqrt;
      using std::isnan;
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;

      if (n == 0)
         return 0;

      const real_t one(1.0);
      const complex<real_t> onec(1.0);
      const real_t zero(0.0);
      real_t ajj;
      complex<real_t> * Aj = A;
      int_t info = 0;
      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t j = 0; j < n; ++j)
         {
            ajj = real(Aj[j] - LATL::dotc(j, Aj, 1, Aj, 1));
            if (ajj <= zero or isnan(ajj))
            {
               Aj[j] = ajj;
               return j+1;
            }
            Aj[j] = ajj = sqrt(ajj);

            if (j < n-1)
            {
               LATL::lacgv(j, Aj, 1);
               LATL::gemv('T', j, n-j-1, -onec, Aj+ldA, ldA, Aj, 1, onec, Aj+ldA+j, ldA);
               LATL::lacgv(j, Aj, 1);
               LATL::scal(n-j-1, one/ajj, Aj+ldA+j, ldA);
            }
            Aj += ldA;
         }
      }
      else
      {
         for (int_t j = 0; j < n; ++j)
         {
            ajj = real(Aj[j] - LATL::dotc(j, A+j, ldA, A+j, ldA));
            if (ajj <= zero || isnan(ajj))
            {
               Aj[j] = ajj;
               return j+1;
            }
            Aj[j] = ajj = sqrt(ajj);

            if (j < n-1)
            {
               LATL::lacgv(j, A+j, ldA);
               LATL::gemv('N', n-j-1, j, -onec, A+j+1, ldA, A+j, ldA, onec, Aj+j+1, 1);
               LATL::lacgv(j, A+j, ldA);
               LATL::scal(n-j-1, one/ajj, Aj+j+1, 1);
            }
            Aj += ldA;
         }
      }
      return info;
   }


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
   /// @param nb Block size.
   
   template< typename real_t>
   int_t potrf(const char uplo, const int_t n, real_t * const A, const int_t ldA, const int_t nb)
   {
      using std::sqrt;
      using std::isnan;
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      if (nb <= 1 || nb >= n)
         return LATL::potrf(uplo, n, A, ldA);
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
               LATL::syrk('U', 'T', jb, j, -one, Aj, ldA, one, Aj+j, ldA);
               info = LATL::potrf('U', jb, Aj+j, ldA);
               if (info != 0)
                  return info+j;
               if (j+jb < n)
               {
                  LATL::gemm('T', 'N', jb, n-j-jb, j, -one, Aj, ldA, Aj+ldA*jb, ldA, one, Aj+ldA*jb+j, ldA);
                  LATL::trsm('L', 'U', 'T', 'N', jb, n-j-jb, one, Aj+j, ldA, Aj+j+ldA*jb, ldA);
               }
               Aj += ldA*nb;
            }
         }
         else
         {
            for (int_t j = 0; j < n; j += nb)
            {
               jb = std::min(nb, n-j);
               LATL::syrk('L', 'N', jb, j, -one, A+j, ldA, one, Aj+j, ldA);
               info = LATL::potrf('L', jb, Aj+j, ldA);
               if (info != 0)
                  return info+j;
               if (j+jb < n)
               {
                  LATL::gemm('N', 'T', n-j-jb, jb, j, -one, A+j+jb, ldA, A+j, ldA, one, Aj+j+jb, ldA );
                  LATL::trsm('R', 'L', 'T', 'N', n-j-jb, jb, one, Aj+j, ldA, Aj+j+jb, ldA);
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
   /// @param nb Block size.

   template< typename real_t>
   int_t potrf(const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA, const int_t nb)
   {
      using std::sqrt;
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      if (nb <= 1 || nb >= n)
         return LATL::potrf(uplo, n, A, ldA);
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
               LATL::herk('U', 'C', jb, j, -one, Aj, ldA, one, Aj+j, ldA);
               info = LATL::potrf('U', jb, Aj+j, ldA);
               if (info != 0)
                  return info+j;
               if (j+jb < n)
               {
                  LATL::gemm('C', 'N', jb, n-j-jb, j, -onec, Aj, ldA, Aj+ldA*jb, ldA, onec, Aj+ldA*jb+j, ldA);
                  LATL::trsm('L', 'U', 'C', 'N', jb, n-j-jb, onec, Aj+j, ldA, Aj+j+ldA*jb, ldA);
               }
               Aj += ldA*nb;
            }
         }
         else
         {
            for (int_t j = 0; j < n; j += nb)
            {
               jb = std::min(nb, n-j);
               LATL::herk('L', 'N', jb, j, -one, A+j, ldA, one, Aj+j, ldA);
               info = LATL::potrf('L', jb, Aj+j, ldA);
               if (info != 0)
                  return info+j;
               if (j+jb < n)
               {
                  LATL::gemm('N', 'C', n-j-jb, jb, j, -onec, A+j+jb, ldA, A+j, ldA, onec, Aj+j+jb, ldA );
                  LATL::trsm('R', 'L', 'C', 'N', n-j-jb, jb, onec, Aj+j, ldA, Aj+j+jb, ldA);
               }
               Aj += ldA*nb;
            }
         }
         return 0;
      }
   }
}

#endif
