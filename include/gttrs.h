//
//  gttrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/22/13.
//
//

#ifndef _gttrs_h
#define _gttrs_h

#include "latl.h"
#include "gtts2.h"

namespace LATL
{
   /// @ingroup COMP
   
   template< typename real_t>
   int_t GTTRS(const char trans, const int_t n, const int_t nrhs, real_t * const DL, real_t * const D, real_t * const DU, real_t * const DU2, int_t * const IPIV, real_t * const B, const int_t ldB, const int_t nb)
   {
      if (trans != 'N' && trans != 'T' && trans != 'C' && trans != 'n' && trans != 't' && trans != 'c')
         return -1;
      if ( n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldB < n)
         return -10;
      if (n == 0 || nrhs == 0)
         return 0;
      int_t itrans = 1;
      using std::min;
      if (trans == 'N' || trans == 'n')
         itrans = 0;
      if (nb >= nrhs)
         LATL::GTTS2(itrans, n, nrhs, DL, D, DU, DU2, IPIV, B, ldB);
      else
      {
         real_t * Bj = B;
         for (int_t j = 0; j < nrhs; j += nb)
         {
            int_t jb = min(nrhs-j, nb);
            LATL::GTTS2(itrans, n, jb, DL, D, DU, DU2, IPIV, Bj, ldB);
            Bj += ldB*nb;
         }
      }
   }
}


#endif
