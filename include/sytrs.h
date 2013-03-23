//
//  sytrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 2/22/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _sytrs_h
#define _sytrs_h

/// @file sytrs.h Solves a system of linear equations A*X = B.

#include "scal.h"
#include "syconv.h"
#include "swap.h"
#include "trsm.h"
#include "latl.h"

namespace LATL
{
   /// @brief Solves a system of linear equations A * X = B with a symmetric n-by-n matrix A using the U*D*U' (if uplo = 'U') or L'*D*L factorization computed by sytrf.
   /// @return 0 
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param A Real matrix size ldA-by-n.  On entry, should contain factors L (or U) and D from the factorization A computed by sytrf.
   /// @param ldA Column length of matrix A.  ldA >= n
   /// @param IPIV Integer array size n.  On entry, contains the details of the interchanges.
   /// @param BSDV Bool array size n.  On entry, contains the details of the block structure of D.
   /// @param B Real matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of B.  ldB >= n
   /// @ingroup SOLV
   
   template<typename real_t>
   int_t sytrs(const char uplo, const int_t n, const int_t nrhs, real_t * const A, const int_t ldA, int_t * IPIV, bool * BSDV, real_t * const B, const int_t ldB)
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
         return -8;
      
      if (n == 0 || nrhs == 0)
         return 0;
      real_t * Work = new real_t[n];
      LATL::syconv(uplo, 'C', n, A, ldA, IPIV, BSDV, Work);
      const real_t one(1.0);
      if (uplo == 'U' || uplo == 'u')
      {
         int_t k = n-1, kp;
         while (k >= 0)
         {
            kp = IPIV[k];
            if (BSDV[k] == 0)
            {
               if (kp != k)
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               --k;
            }
            else
            {
               if (kp == IPIV[k-1])
               {
                  LATL::swap(nrhs, B+k-1, ldB, B+kp, ldB);
               }
               k -= 2;
            }
         }
         LATL::trsm('L', 'U', 'N', 'U', n, nrhs, one, A, ldA, B, ldB);
         
         int_t i = n-1;
         real_t * Ai;
         while (i >= 0)
         {
            Ai = A+ldA*i;
            if (BSDV[i] == 0)
            {
               LATL::scal(nrhs, one/Ai[i], B+i, ldB);
            }
            else if (i > 0)
            {
               if (IPIV[i-1] == IPIV[i])
               {
                  real_t * Aim1 = Ai - ldA;
                  real_t akm1k = Work[i];
                  real_t akm1 = Aim1[i-1]/akm1k;
                  real_t ak = Ai[i]/akm1k;
                  real_t denom = akm1*ak-one;
                  real_t * Bj = B;
                  for (int_t j = 0; j < nrhs; ++j)
                  {
                     real_t bkm1 = Bj[i-1]/akm1k;
                     real_t bk = Bj[i]/akm1k;
                     Bj[i-1] = (ak*bkm1-bk)/denom;
                     Bj[i] = (akm1*bk-bkm1)/denom;
                     Bj += ldB;
                  }
                  --i;
               }
            }
            --i;
         }
         trsm('L', 'U', 'T', 'U', n, nrhs,  one, A, ldA, B, ldB);
         k = 0;
         while (k < n)
         {
            kp = IPIV[k];
            if (BSDV[k] == 0)
            {
               if (kp != k)
               {
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               ++k;
            }
            else
            {
               if (k < n-1 && kp == IPIV[k+1])
               {
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               k+=2;
            }
         }
      }
      else
      {
         int_t k = 0, kp;
         while (k < n)
         {
            if (BSDV[k] == 0)
            {
               kp = IPIV[k];
               if (kp != k)
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               ++k;
            }
            else
            {
               kp = IPIV[k+1];
               if (kp == IPIV[k])
               {
                  LATL::swap(nrhs, B+k+1, ldB, B+kp, ldB);
               }
               k += 2;
            }
         }
         LATL::trsm('L', 'L', 'N', 'U', n, nrhs, one, A, ldA, B, ldB);
         int_t i = 0;
         while (i < n)
         {
            real_t * Ai = A+ldA*i;
            if (BSDV[i] == 0)
            {
               LATL::scal(nrhs, one/Ai[i], B+i, ldB);
            }
            else
            {
               real_t * Aip1 = Ai+ldA;
               real_t akm1k = Work[i];
               real_t akm1 = Ai[i]/akm1k;
               real_t ak = Aip1[i+1]/akm1k;
               real_t denom = akm1*ak-one;
               real_t * Bj = B;
               for (int_t j = 0; j < nrhs; ++j)
               {
                  real_t bkm1 = Bj[i]/akm1k;
                  real_t bk = Bj[i+1]/akm1k;
                  Bj[i] = (ak*bkm1-bk)/denom;
                  Bj[i+1] = (akm1*bk-bkm1)/denom;
                  Bj += ldB;
               }
               ++i;
            }
            ++i;
         }
         LATL::trsm('L', 'L', 'T', 'U', n, nrhs, one, A, ldA, B, ldB);
         k = n-1;
         while (k >= 0)
         {
            kp = IPIV[k];
            if (BSDV[k] == 0)
            {
               if (kp != k)
               {
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               --k;
            }
            else
            {
               if (k > 0 && kp == IPIV[k-1])
               {
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               k -=2;
            }
         }
      }
   
      LATL::syconv(uplo, 'R', n, A, ldA, IPIV, BSDV, Work);
      delete [] Work;
      return 0;
   }

   /// @brief Solves a system of linear equations A * X = B with a complex symmetric n-by-n matrix A using the U*D*U' (if uplo = 'U') or L'*D*L factorization computed by sytrf.
   /// @return 0
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param A Complex matrix size ldA-by-n.  On entry, should contain factors L (or U) and D from the factorization A computed by sytrf.
   /// @param ldA Column length of matrix A.  ldA >= n
   /// @param IPIV Integer array size n.  On entry, contains the details of the interchanges.
   /// @param BSDV Bool array size n.  On entry, contains the details of the block structure of D.
   /// @param B Real matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of B.  ldB >= n
   /// @ingroup SOLV
   
   template<typename real_t>
   int_t sytrs(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, const int_t ldA, int_t * const IPIV, bool * const BSDV, complex<real_t> * const B, const int_t ldB)
   {
      return LATL::sytrs< complex<real_t> >(uplo, n, nrhs, A, ldA, IPIV, BSDV, B, ldB);
   }
}

#endif
