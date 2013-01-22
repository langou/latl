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
   ///
   /// @tparam real_t Floating point type.
   /// @param n Number of columns in matrix A.
   /// @param A Real matrix size ldA-by-n.  On exit, the permuted matrix.
   /// @param ldA Column length of A.
   /// @param k1 The index of the first element of IPIV for a row interchange.
   /// @param k2 The index of the final element of IPIV for a row interchange.
   /// @param IPIV Integer array length at least k2+1.  For each k between k1 and k2, IPIV[k] = L indicates an exchange of row k of A with row L.
   /// @param forward Boolean value, optional.  Defaults to 1 for reading indices from the first index, k1, to the last, k2.  A value of 0 will read the values of IPIV in reverse order, from k2 to k1.
   /// @ingroup MAT
   
   template< typename real_t>
   void laswp(const int_t n, real_t * const A, const int_t ldA, const int_t k1, const int_t k2, int_t * const IPIV, const bool forward = 1)
   {
      int_t i1, i2, inc;
      if (forward == 1)
      {
         i1 = k1;
         i2 = k2;
         inc = 1;
      }
      else
      {
         i1 = k2;
         i2 = k1;
         inc = -1;
      }
      
      int_t n32 = (n/32)*32;
      int_t ix, ip, colSwap;
      real_t temp;
      real_t *Ai = A + i1, *Aip = A;
      if (n32 != 0)
      {
         for (int_t j = 0; j < n32; j+= 32)
         {
            ix = i1;
            for (int_t i = i1; i <= i2; i += inc)
            {
               ip = IPIV[ix];
               if (ip != i)
               {
                  Aip = A + ip;
                  for (int_t k = j; k < j+32; ++k)
                  {
                     int_t colSwap = ldA*k;
                     temp = Ai[colSwap];
                     Ai[colSwap] = Aip[colSwap];
                     Aip[colSwap] = temp;
                  }
               }
               ix += inc;
               Ai += inc;
            }
         }
      }
      if (n32 != n)
      {
         ix = i1;
         Ai = A;
         for (int_t i = i1; i <= i2; i+= inc)
         {
            ip = IPIV[ix];
            if (ip != i)
            {
               Aip = A + ip;
               for (int_t k = n32; k < n; ++k)
               {
                  colSwap = ldA*k;
                  temp = Ai[colSwap];
                  Ai[colSwap] = Aip[colSwap];
                  Aip[colSwap] = temp;
               }
            }
            ix += inc;
            Ai += inc;
         }
      }
      return;
   }
   
   /// @brief Applies a series of row changes to the matrix A.
   ///
   /// @tparam real_t Floating point type.
   /// @param n Number of columns in matrix A.
   /// @param A Complex matrix size ldA-by-n.  On exit, the permuted matrix.
   /// @param ldA Column length of A.
   /// @param k1 The index of the first element of IPIV for a row interchange.
   /// @param k2 The index of the final element of IPIV for a row interchange.
   /// @param IPIV Integer array length at least k2+1.  For each k between k1 and k2, IPIV[k] = L indicates an exchange of row k of A with row L.
   /// @param forward Boolean value, optional.  Defaults to 1 for reading indices from the first index, k1, to the last, k2.  A value of 0 will read the values of IPIV in reverse order, from k2 to k1.
   /// @ingroup MAT
   
   template< typename real_t>
   void laswp(const int_t n, complex<real_t> * const A, const int_t ldA, const int_t k1, const int_t k2, int_t * const IPIV, const bool forward = 1)
   {
      int_t i1, i2, inc;
      if (forward == 1)
      {
         i1 = k1;
         i2 = k2;
         inc = 1;
      }
      else
      {
         i1 = k2;
         i2 = k1;
         inc = -1;
      }
      
      int_t n32 = (n/32)*32;
      int_t ix, ip, colSwap;
      complex<real_t> temp;
      complex<real_t> *Ai = A + i1, *Aip = A;
      if (n32 != 0)
      {
         for (int_t j = 0; j < n32; j+= 32)
         {
            ix = i1;
            for (int_t i = i1; i <= i2; i += inc)
            {
               ip = IPIV[ix];
               if (ip != i)
               {
                  Aip = A + ip;
                  for (int_t k = j; k < j+32; ++k)
                  {
                     int_t colSwap = ldA*k;
                     temp = Ai[colSwap];
                     Ai[colSwap] = Aip[colSwap];
                     Aip[colSwap] = temp;
                  }
               }
               ix += inc;
               Ai += inc;
            }
         }
      }
      if (n32 != n)
      {
         ix = i1;
         Ai = A;
         for (int_t i = i1; i <= i2; i+= inc)
         {
            ip = IPIV[ix];
            if (ip != i)
            {
               Aip = A + ip;
               for (int_t k = n32; k < n; ++k)
               {
                  colSwap = ldA*k;
                  temp = Ai[colSwap];
                  Ai[colSwap] = Aip[colSwap];
                  Aip[colSwap] = temp;
               }
            }
            ix += inc;
            Ai += inc;
         }
      }
      return;
   }
}
#endif
