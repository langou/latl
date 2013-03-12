//
//  sysv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 2/19/13.
//
//

#ifndef _sysv_h
#define _sysv_h

/// @file sysv.h

#include "sytrf.h"
#include "sytrs.h"

namespace latl
{
   template<typename real_t>
   int_t sysv(const char uplo, const int_t n, const int_t nrhs, real_t * const A, int_t ldA, int_t * ipiv, bool * bsdv, real_t * const B, int_t ldB)
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
      
      int_t info = latl::sytrf(uplo, n, A, ldA, ipiv, bsdv);
      if (info == 0)
      {
         info = latl::sytrs(uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
      }
      return info;
   }
   
   template<typename real_t>
   int_t sysv(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, int_t ldA, int_t * ipiv, bool * bsdv, complex<real_t> * const B, int_t ldB)
   {
      return latl::sysv< complex<real_t> > (uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
   }
   
   template<typename real_t>
   int_t sysv(const char uplo, const int_t n, const int_t nrhs, real_t * const A, int_t ldA, int_t * ipiv, bool * bsdv, real_t * const B, int_t ldB, int_t nb = 32)
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
      
      int_t info = latl::sytrf(uplo, n, A, ldA, ipiv, bsdv, nb);
      if (info == 0)
      {
         info = latl::sytrs(uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
      }
      return info;
   }
   
   template<typename real_t>
   int_t sysv(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, int_t ldA, int_t * ipiv, bool * bsdv, complex<real_t> * const B, int_t ldB, int_t nb = 32)
   {
      return latl::sysv< complex<real_t> > (uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB, nb);
   }
}
#endif
