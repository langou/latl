//
//  hetrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/11/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _hetrs_h
#define _hetrs_h

/// @file hetrs.h Solves a system of linear equations A*X = B.

#include "scal.h"
#include "syconv.h"
#include "swap.h"
#include "trsm.h"
#include "latl.h"

namespace LATL
{
   /// @brief Solves a system of linear equations A * X = B with a Hermitian n-by-n matrix A using the U*D*U^H (if uplo = 'U') or L^H*D*L factorization computed by sytrf.
   /// @return 0
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the Hermitian matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param A Complex Hermitian matrix size ldA-by-n.  On entry, should contain factors L (or U) and D from the factorization A computed by hetrf.
   /// @param ldA Column length of matrix A.  ldA >= n
   /// @param IPIV Integer array size n.  On entry, contains the details of the interchanges.
   /// @param BSDV Bool array size n.  On entry, contains the details of the block structure of D.
   /// @param B Real matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of B.  ldB >= n
   /// @ingroup SOLV
   
   template<typename real_t>
   int_t hetrs(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, const int_t ldA, int_t * ipiv, bool * bsdv, complex<real_t> * const B, const int_t ldB)
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
      complex<real_t> * Work = new complex<real_t>[n];
      LATL::syconv(uplo, 'C', n, A, ldA, ipiv, bsdv, Work);
      const complex<real_t> one(1.0);
      if (uplo == 'U' || uplo == 'u')
      {
         int_t k = n-1, kp;
         while (k >= 0)
         {
            kp = ipiv[k];
            if (bsdv[k] == 0)
            {
               if (kp != k)
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               --k;
            }
            else
            {
               if (kp == ipiv[k-1])
               {
                  LATL::swap(nrhs, B+k-1, ldB, B+kp, ldB);
               }
               k -= 2;
            }
         }

         LATL::trsm('L', 'U', 'N', 'U', n, nrhs, one, A, ldA, B, ldB);
         int_t i = n-1;
         complex<real_t> * Ai;
         while (i >= 0)
         {
            Ai = A+ldA*i;
            if (bsdv[i] == 0)
            {
               real_t s = real(one)/real(Ai[i]);
               LATL::scal(nrhs, s, B+i, ldB);
            }
            else if (i > 0)
            {
               if (ipiv[i-1] == ipiv[i])
               {
                  complex<real_t> * Aim1 = Ai - ldA;
                  complex<real_t> akm1k = Work[i];
                  complex<real_t>akm1 = Aim1[i-1]/akm1k;
                  complex<real_t>ak = Ai[i]/conj(akm1k);
                  complex<real_t>denom = akm1*ak-one;
                  complex<real_t>* Bj = B;
                  for (int_t j = 0; j < nrhs; ++j)
                  {
                     complex<real_t>bkm1 = Bj[i-1]/akm1k;
                     complex<real_t>bk = Bj[i]/conj(akm1k);
                     Bj[i-1] = (ak*bkm1-bk)/denom;
                     Bj[i] = (akm1*bk-bkm1)/denom;
                     Bj += ldB;
                  }
                  --i;
               }
            }
            --i;
         }
         trsm('L', 'U', 'C', 'U', n, nrhs,  one, A, ldA, B, ldB);
         k = 0;
         while (k < n)
         {
            kp = ipiv[k];
            if (bsdv[k] == 0)
            {
               if (kp != k)
               {
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               ++k;
            }
            else
            {
               if (k < n-1 && kp == ipiv[k+1])
               {
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               k+=2;
            }
         }
      }
      else // lower
      {
         int_t k = 0, kp;
         while (k < n)
         {
            if (bsdv[k] == 0)
            {
               kp = ipiv[k];
               if (kp != k)
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               ++k;
            }
            else
            {
               kp = ipiv[k+1];
               if (kp == ipiv[k] && kp != k+1)
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
            complex<real_t>* Ai = A+ldA*i;
            if (bsdv[i] == 0)
            {
               real_t s = real(one)/real(Ai[i]);
               LATL::scal(nrhs, s, B+i, ldB);
               ++i;
            }
            else
            {
               complex<real_t>* Aip1 = Ai+ldA;
               complex<real_t> akm1k = Work[i];
               complex<real_t> akm1 = Ai[i]/conj(akm1k);
               complex<real_t> ak = Aip1[i+1]/akm1k;
               complex<real_t> denom = akm1*ak-one;
               complex<real_t>* Bj = B;
               for (int_t j = 0; j < nrhs; ++j)
               {
                  complex<real_t> bkm1 = Bj[i]/conj(akm1k);
                  complex<real_t> bk = Bj[i+1]/akm1k;
                  Bj[i] = (ak*bkm1-bk)/denom;
                  Bj[i+1] = (akm1*bk-bkm1)/denom;
                  Bj += ldB;
               }
               i +=2;
            }
         }
         LATL::trsm('L', 'L', 'C', 'U', n, nrhs, one, A, ldA, B, ldB);
         k = n-1;
         while (k >= 0)
         {
            kp = ipiv[k];
            if (bsdv[k] == 0)
            {
               if (kp != k)
               {
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               --k;
            }
            else
            {
               if (k > 0 && kp == ipiv[k-1])
               {
                  LATL::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               k -=2;
            }
         }
      }
      
      LATL::syconv(uplo, 'R', n, A, ldA, ipiv, bsdv, Work);
      delete [] Work;
      return 0;
   }
}


#endif
