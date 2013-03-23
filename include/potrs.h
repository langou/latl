//
//  potrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/12/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _potrs_h
#define _potrs_h

/// @file potrs.h

#include "trsm.h"
#include "latl.h"

namespace LATL
{
   template< typename real_t>
   int_t POTRS( const char uplo, const int_t n, const int_t nrhs, real_t * const A, const int_t ldA, real_t * const B, const int_t ldB)
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
      
      if (n == 0 || nrhs == 0)
         return 0;
      const real_t one(1.0);
      if (uplo == 'U' || uplo == 'u')
      {
         LATL::TRSM('L', 'U', 'T', 'N', n, nrhs, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'U', 'N', 'N', n, nrhs, one, A, ldA, B, ldB);
      }
      else
      {
         LATL::TRSM('L', 'L', 'N', 'N', n, nrhs, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'L', 'T', 'N', n, nrhs, one, A, ldA, B, ldB);
      }
      return 0;
   }
   
   template< typename real_t>
   int_t POTRS( const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, const int_t ldA, complex<real_t> * const B, const int_t ldB)
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
      
      if (n == 0 || nrhs == 0)
         return 0;
      const complex<real_t> one(1.0);
      if (uplo == 'U' || uplo == 'u')
      {
         LATL::TRSM('L', 'U', 'C', 'N', n, nrhs, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'U', 'N', 'N', n, nrhs, one, A, ldA, B, ldB);
      }
      else
      {
         LATL::TRSM('L', 'L', 'N', 'N', n, nrhs, one, A, ldA, B, ldB);
         LATL::TRSM('L', 'L', 'C', 'N', n, nrhs, one, A, ldA, B, ldB);
      }
      return 0;
   }
}

#endif
