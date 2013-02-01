//
//  laswp.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 5/16/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _laswp_h
#define _laswp_h

/// @file laswp.h Applies a series of row changes to the matrix A.

#include "latl.h"

namespace latl
{
   /// @brief Applies a series of row changes to the matrix A.
   /// @tparam real_t Floating point type.
   /// @param n Number of columns in matrix A.
   /// @param A Real matrix size ldA-by-n.  On exit, the permuted matrix.
   /// @param ldA Column length of A.
   /// @param k1 The index of the first element of IPIV for a row interchange.
   /// @param k2 The index of the final element of IPIV for a row interchange.
   /// @param IPIV Integer array length at least k2+1.  For each k between k1 and k2,
   /// IPIV[k] = L indicates an exchange of row k of A with row L.
   /// @param inc Determines whether IPIV is read forward (inc=1) or backward (inc=-1).
   /// (optional, default value is 1)
   /// @ingroup MAT
   
   template< typename real_t>
   int laswp(const int_t n, real_t * const A, const int_t ldA, const int_t k1, const int_t k2, int_t * const IPIV, int_t inc=1)
   {
      const int_t b=32;
      int_t i1, i2, i0;
      if (inc==1)
      {
         i0 = k1;
         i1 = k1;
         i2 = k2;
      }
      else if (inc==-1)
      {
         i0 = k2;
         i1 = k2;
         i2 = k1;
      }
      else
         return 0;
      
      int_t nb = (n/b)*b;
      real_t *Aj = A;
      int_t ip,ix;
      for (int_t j = 0; j < nb; j+= b)
      {
         ix = i0;
         for (int_t i = i1; i <= i2; i += inc)
         {
            ip = IPIV[ix];
            if (ip != i)
            {
               real_t *Ak = Aj;
               for (int_t k = j; k < j+b; k++)
               {
                  real_t temp = Ak[i];
                  Ak[i] = Ak[ip];
                  Ak[ip] = temp;
                  Ak+=ldA;
               }
            }
            ix += inc;
         }
         Aj += b*ldA;
      }
      if (nb != n)
      {
         ix = i0;
         for (int_t i = i1; i <= i2; i+= inc)
         {
            ip = IPIV[ix];
            if (ip != i)
            {
               real_t *Ak = A+nb*ldA;
               for (int_t k = nb; k < n; k++)
               {
                  real_t temp = Ak[i];
                  Ak[i] = Ak[ip];
                  Ak[ip] = temp;
                  Ak+=ldA;
               }
            }
            ix += inc;
         }
      }
      return 0;
   }
   
   /// @brief Applies a series of row changes to the matrix A.
   /// @tparam real_t Floating point type.
   /// @param n Number of columns in matrix A.
   /// @param A Complex matrix size ldA-by-n.  On exit, the permuted matrix.
   /// @param ldA Column length of A.
   /// @param k1 The index of the first element of IPIV for a row interchange.
   /// @param k2 The index of the final element of IPIV for a row interchange.
   /// @param IPIV Integer array length at least k2+1.  For each k between k1 and k2,
   /// IPIV[k] = L indicates an exchange of row k of A with row L.
   /// @param inc Determines whether IPIV is read forward (inc=1) or backward (inc=-1).
   /// (optional, default value is 1)
   /// @ingroup MAT
   
   template< typename real_t>
   int laswp(const int_t n, complex<real_t> * const A, const int_t ldA, const int_t k1, const int_t k2, int_t * const IPIV, int_t inc=1)
   {
      return laswp< complex<real_t> >(n,A,ldA,k1,k2,IPIV,inc);
   }
}
#endif
