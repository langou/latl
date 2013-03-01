//
//  sytrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 2/19/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _sytrs_h
#define _sytrs_h

/// @file sytrs.h

#include "gemv.h"
#include "ger.h"
#include "scal.h"
#include "swap.h"
#include "latl.h"

namespace latl
{
   template<typename real_t>
   int_t sytrs(const char uplo, const int_t n, const int_t nrhs, real_t * const A, const int_t ldA, int_t * const IPIV, bool * const BSDV, real_t * const B, const int_t ldB)
   {
      if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
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
      const real_t one(1.0);
      if (uplo == 'U' || uplo == 'u')
      {
         int_t k = n-1, kp;
         real_t *Ak;
         while (k >=0)
         {
            Ak = A + k*ldA;
            if (BSDV[k] == 0)
            {
               kp = IPIV[k];
               if (kp != k)
               {
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               latl::ger(k, nrhs, -one, A+k, 1, B+k, ldB, B, ldB);
               latl::scal(nrhs, one/Ak[k], B+k, ldB);
               --k;
            }
            else
            {
               real_t * Akm1 = Ak-ldA;
               kp = IPIV[k];
               if (kp != k-1)
               {
                  latl::swap(nrhs, B+k-1, ldB, B+kp, ldB);
               }
               latl::ger(k-1, nrhs, -one, Ak, 1, B+k, ldB, B, ldB);
               latl::ger(k-1, nrhs, -one, Akm1, 1, B+k-1, ldB, B, ldB);
               
               real_t akkm1 = Ak[k-1];
               real_t akm1 = Akm1[k-1]/akkm1;
               real_t ak = Ak[k]/akkm1;
               real_t denom = akm1*ak-one;
               real_t * Bj = B;
               for (int_t j = 0; j < nrhs; ++j)
               {
                  
                  real_t bkm1 = Bj[k-1]/akkm1;
                  real_t bk = Bj[k]/akkm1;
                  Bj[k-1] = (ak*bkm1-bk)/denom;
                  Bj[k] = (akm1*bk-bkm1)/denom;
                  Bj += ldB;
               }
            k -= 2;
            }
         }

         k=0;
         while (k < n)
         {
            Ak = A+ldA*k;
            if (BSDV[k] == 0)
            {
               latl::gemv('T', k, nrhs, -one, B, ldB, Ak, 1, one, B+k, ldB);
               kp = IPIV[k];
               if (kp != k)
               {
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               ++k;
            }
            else
            {
               latl::gemv('T', k, nrhs, -one, B, ldB, Ak, 1, one, B+k, ldB);
               latl::gemv('T', k, nrhs, -one, B, ldB, Ak+ldA, 1, one, B+k+1, ldB);
               kp = IPIV[k];
               if (kp != k)
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               k += 2;
            }
         }
         
      }
      else
      {
         int_t k = 0, kp;
         real_t *Ak;
         while (k < n)
         {
            Ak = A+ldA*k;
            if (BSDV[k] == 0)
            {
               kp = IPIV[k];
               if (kp != k)
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               if (k < n-1)
                  latl::ger(n-k-1, nrhs, -one, Ak+k+1, 1, B+k, ldB, B+k+1, ldB);
               latl::scal(nrhs, one/Ak[k], B+k, ldB);
               ++k;
            }
            else
            {
               real_t * Akp1 = Ak+ldA;
               kp = IPIV[k];
               if (kp != k+1)
                  latl::swap(nrhs, B+k+1, ldB, B+kp, ldB);
               if (k < n-2)
               {
                  latl::ger(n-k-2, nrhs, -one, Ak+k+2, 1, B+k, ldB, B+k+2, ldB);
                  latl::ger(n-k-2, nrhs, -one, Akp1+k+2, 1, B+k+1, ldB, B+k+2, ldB);
               }
               
               real_t akkp1 = Ak[k+1];
               real_t akk = Ak[k]/akkp1;
               real_t ak = Akp1[k+1]/akkp1;
               real_t denom = akk*ak-one;
               real_t * Bj = B;
               for (int_t j = 0; j < nrhs; ++j)
               {
                  real_t bk = Bj[k]/akkp1;
                  real_t bkp1 = Bj[k+1]/akkp1;
                  Bj[k] = (ak*bk-bkp1)/denom;
                  Bj[k+1] = (akk*bkp1-bk)/denom;
                  Bj += ldB;
               }
               k += 2;
            }
         }
         k = n-1;
         while (k >= 0)
         {
            Ak = A+ldA*k;
            if (BSDV[k] == 0)
            {
               if (k < n-1)
               {
                  latl::gemv('T', n-k-1, nrhs, -one, B+k+1, ldB, Ak+k+1, 1, one, B+k, ldB);
               }
               kp = IPIV[k];
               if (kp != k)
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               --k;
            }
            else
            {
               if (k < n-1)
               {
                  latl::gemv('T', n-k-1, nrhs, -one, B+k+1, ldB, Ak+k+1, 1, one, B+k, ldB);
                  latl::gemv('T', n-k-1, nrhs, -one, B+k+1, ldB, Ak-ldA+k+1, 1, one, B+k-1, ldB);
               }
               kp = IPIV[k];
               if (kp != k)
               {
                  latl::swap(nrhs, B+k, ldB, B+kp, ldB);
               }
               k -= 2;
            }
         }
      }
      return 0;
   }
   
   template<typename real_t>
   int_t sytrs(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, const int_t ldA, int_t * const IPIV, bool * const BSDV, complex<real_t> * const B, const int_t ldB)
   {
      return latl::sytrs< complex<real_t> >(uplo, n, nrhs, A, ldA, IPIV, BSDV, B, ldB);
   }
}


#endif
