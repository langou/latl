//
//  pttrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/21/13.
//
//

#ifndef _pttrs_h
#define _pttrs_h

/// @file pttrs.h

#include "latl.h"
#include "ptts2.h"

namespace latl
{
   template< typename real_t>
   int_t pttrs(const int_t n, const int_t nrhs, real_t * const D, real_t * const E, real_t * const B, const int_t ldB, const int_t nb)
   {
      if ( n < 0)
         return -1;
      if ( nrhs < 0)
         return -2;
      if (ldB < n)
         return -6;
      
      if ( n == 0 || nrhs == 0)
         return 0;
      using std::min;
      if (nb >= nrhs)
         latl::ptts2(n, nrhs, D, E, B, ldB);
      else
      {
         real_t * Bj = B;
         for (int_t j = 0; j < nrhs; j += nb)
         {
            int_t jb = min(nrhs-j, nb);
            latl::ptts2(n, jb, D, E, Bj, ldB);
            Bj += ldB*nb;
         }
      }
   }
}

#endif
