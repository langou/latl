//
//  ptsv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/21/13.
//
//

#ifndef _ptsv_h
#define _ptsv_h

/// @file ptsv.h

#include "latl.h"
#include "pttrs.h"

namespace latl
{
   template< typename real_t>
   int_t ptsv(const int_t n, const int_t nrhs, real_t * const D, real_t * const E, real_t * const B, const int_t ldB, const int_t nb)
   {
      if ( n < 0)
         return -1;
      if (nrhs < 0)
         return -2;
      if (ldB<n)
         return -6;
      
      int_t info;
      info = latl::pttrf(n, D, E);
      if (info == 0)
      {
         info = latl::pttrs(n, nrhs, D, E, B, ldB, nb);
      }
      return info;
   }
}


#endif
