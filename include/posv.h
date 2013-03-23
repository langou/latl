//
//  posv.h
//  
//
//  Created by Stephanie Patterson on 3/14/13.
//
//

#ifndef _posv_h
#define _posv_h

#include "potrs.h"
#include "potrf.h"

namespace LATL
{
   template<typename real_t>
   int_t posv(const char uplo, const int_t n, const int_t nrhs, real_t * const A, const int_t ldA, real_t * const B, const int_t ldB)
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
         return -7;
      int_t info = 0;
      info = LATL::potrf(uplo, n, A, ldA);
      if (info == 0)
          info = LATL::potrs(uplo, n, nrhs, A, ldA, B, ldB);
      return info;
   }
   
   template<typename real_t>
   int_t posv(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, const int_t ldA, complex<real_t> * const B, const int_t ldB)
   {
      return LATL::posv< complex<real_t> >(uplo, n, nrhs, A, ldA, B, ldB);
   }
}

#endif
