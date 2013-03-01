//
//  sytrs2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 2/22/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _sytrs2_h
#define _sytrs2_h

/// @file sytrs2.h

#include "scal.h"
#include "syconv.h"
#include "swap.h"
#include "trsm.h"
#include "latl.h"

namespace latl
{
   template<typename real_t>
   int_t sytrs2(const char uplo, const int_t n, const int_t nrhs, real_t * const A, const int_t ldA, int_t * ipiv, bool * bsdv, real_t * const B, const int_t ldB)
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
      latl::syconv(uplo, 'C', n, A, ldA, ipiv, bsdv, Work);
      const real_t one(1.0);
      if (uplo == 'U' || uplo == 'u')
      {
         int_t k = n-1, kp;
         while (k >= 0)
         {
            kp = ipiv[k];
            if (bsdv[k] == 0)
            {
               if (kp != k)
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               --k;
            }
            else
            {
               if (kp == ipiv[k-1])
               {
                  latl::swap(nrhs, B+k-1, ldB, B+kp, ldB);
               }
               k -= 2;
            }
         }
         latl::trsm('L', 'U', 'N', 'U', n, nrhs, one, A, ldA, B, ldB);
         
         int_t i = n-1;
         real_t * Ai;
         while (i >= 0)
         {
            Ai = A+ldA*i;
            if (bsdv[i] == 0)
            {
               latl::scal(nrhs, one/Ai[i], B+i, ldB);
            }
            else if (i > 0)
            {
               if (ipiv[i-1] == ipiv[i])
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
            kp = ipiv[k];
            if (bsdv[k] == 0)
            {
               if (kp != k)
               {
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               ++k;
            }
            else
            {
               if (k < n-1 && kp == ipiv[k+1])
               {
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
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
            if (bsdv[k] == 0)
            {
               kp = ipiv[k];
               if (kp != k)
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               ++k;
            }
            else
            {
               kp = ipiv[k+1];
               if (kp == ipiv[k])
               {
                  latl::swap(nrhs, B+k+1, ldB, B+kp, ldB);
               }
               k += 2;
            }
         }
         latl::trsm('L', 'L', 'N', 'U', n, nrhs, one, A, ldA, B, ldB);
         int_t i = 0;
         while (i < n)
         {
            real_t * Ai = A+ldA*i;
            if (bsdv[i] == 0)
            {
               latl::scal(nrhs, one/Ai[i], B+i, ldB);
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
         latl::trsm('L', 'L', 'T', 'U', n, nrhs, one, A, ldA, B, ldB);
         k = n-1;
         while (k >= 0)
         {
            kp = ipiv[k];
            if (bsdv[k] == 0)
            {
               if (kp != k)
               {
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               --k;
            }
            else
            {
               if (k > 0 && kp == ipiv[k-1])
               {
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               k -=2;
            }
         }
      }
   
      latl::syconv(uplo, 'R', n, A, ldA, ipiv, bsdv, Work);
      delete [] Work;
      return 0;
   }

   template<typename real_t>
   int_t sytrs2(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, const int_t ldA, int_t * const ipiv, bool * const bsdv, complex<real_t> * const B, const int_t ldB)
   {
      return latl::sytrs2< complex<real_t> >(uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
   }
}

#endif
