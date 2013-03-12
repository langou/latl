//
//  hetrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/11/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _hetrs_h
#define _hetrs_h

/// @file hetrs.h

#include "scal.h"
#include "syconv.h"
#include "swap.h"
#include "trsm.h"
#include "latl.h"

namespace latl
{

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
      latl::syconv(uplo, 'C', n, A, ldA, ipiv, bsdv, Work);
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
         complex<real_t> * Ai;
         while (i >= 0)
         {
            Ai = A+ldA*i;
            if (bsdv[i] == 0)
            {
               real_t s = real(one)/real(Ai[i]);
               latl::scal(nrhs, s, B+i, ldB);
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
      else // lower
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
               if (kp == ipiv[k] && kp != k+1)
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
            complex<real_t>* Ai = A+ldA*i;
            if (bsdv[i] == 0)
            {
               real_t s = real(one)/real(Ai[i]);
               latl::scal(nrhs, s, B+i, ldB);
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
         latl::trsm('L', 'L', 'C', 'U', n, nrhs, one, A, ldA, B, ldB);
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
}


#endif
