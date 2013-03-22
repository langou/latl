//
//  ptts2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/20/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _ptts2_h
#define _ptts2_h

///@file ptts2.h

#include "scal.h"
#include "latl.h"

namespace latl
{
   template<typename real_t>
   int_t ptts2(const int_t n, const int_t nrhs, real_t * const D, real_t * const E, real_t * const B, const int_t ldB)
   {
      if (n < 0)
         return -1;
      if (nrhs < 0)
         return -2;
      if (ldB < n)
         return -6;
      if (n == 0 || nrhs == 0)
         return 0;
      if (n == 1)
      {
         latl::scal(nrhs, 1.0/D[0], B, ldB);
         return 0;
      }
      real_t * Bj = B;
      for (int_t j = 0; j < nrhs; ++j)
      {
         for (int_t i = 1; i < n; ++i)
         {
            Bj[i] -= Bj[i-1]*E[i-1];
         }
         Bj[n-1] /= D[n-1];
         for (int_t i = n-2; i >= 0; --i)
         {
            Bj[i] = Bj[i]/D[i] - Bj[i+1]*E[i];
         }
         Bj += ldB;
      }
   }
   
   template<typename real_t>
   int_t ptts2(bool iuplo, const int_t n, const int_t nrhs, real_t * const D, complex<real_t> * const E, complex<real_t> * const B, const int_t ldB)
   {
      if ( n < 0)
         return -1;
      if (nrhs < 0)
         return -2;
      if (ldB < n)
         return -6;
      if (n == 0 || nrhs == 0)
         return 0;
      if (n == 1)
      {
         latl::scal(nrhs, 1.0/D[0], B, ldB);
         return 0;
      }
      if (iuplo == 1)
      {
         complex<real_t> * Bj = B;
         for (int_t j = 0; j < nrhs; ++j)
         {
            for (int_t i = 1; i < n; ++i)
            {
               Bj[i] -= Bj[i-1]*conj(E[i-1]);
            }
            Bj[n-1] = Bj[n-1]/D[n-1];
            for (int_t i = n-2; i >= 0; --i)
            {
               Bj[i] = Bj[i]/D[i]-Bj[i+1]*E[i];
            }
            Bj += ldB;
         }
      }
      else
      {
         complex<real_t> * Bj = B;
         for (int_t j = 0; j < nrhs; ++j)
         {
            for (int_t i = 1; i < n; ++i)
            {
               Bj[i] -= Bj[i-1]*E[i-1];
            }
            Bj[n-1] = Bj[n-1]/D[n-1];
            for (int_t i = n-2; i >= 0; --i)
            {
               Bj[i] = Bj[i]/D[i]-Bj[i+1]*conj(E[i]);
            }
            Bj += ldB;
         }
      }
   }
}

#endif
