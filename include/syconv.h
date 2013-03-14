//
//  syconv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 2/22/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _syconv_h
#define _syconv_h

/// @file syconv.h Converts A, given by sytrf or hetrf, into L and D and vice-versa.

#include "latl.h"

namespace latl
{
   /// @brief Converts A, given by sytrf, into L and D and vice-versa.
   ///
   /// The off-diagonal elements of D are returned in Work.  Also applies (or reverses) permutation done in TRF.
   ///
   /// This is an auxiliary routine for sytrs.
   /// @return 0
   /// @tparam real_t Floating point type
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.
   /// @param way Indicates whether the matrix A should be converted to L and D, or if the changes need to be reverted.
   ///
   ///      'C' = convert
   ///      'R' = revert
   ///
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Real symmetric matrix size ldA-by-n.  On entry, the block diagonal matrix D and the multipliers used to obtain the factor U or L as computed by sytrf.
   /// @param ldA Column length or the matrix A.  ldA >=n
   /// @param ipiv Integer array size n.  On entry, should contain the details of the interchanges of D.
   /// @param bsdv Bool array size n.  On entry, should contain the details of the block structure of D, as computed by sytrf.
   /// @param Work Real array length n.  On exit, contains the off-diagonal elements of A, with any other elements set to zero.
   /// @ingroup AUX
   
   template< typename real_t>
   int_t syconv(const char uplo, const char way, const int_t n, real_t * const A, const int_t ldA, int_t * ipiv, bool * bsdv, real_t * Work)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (way != 'C' && way != 'R' && way != 'c' && way != 'r')
         return -2;
      if ( n < 0)
         return -3;
      if (ldA < n)
         return -5;
      
      if (n == 0)
         return 0;
      
      int_t i = 0, ip;
      const real_t zero(0.0);
      if (uplo == 'U' || uplo == 'u')
      {
         if (way == 'C' || way == 'c')
         {
            i = n-1;
            Work[0] = zero;
            real_t * Ai;
            while (i > 0)
            {
               if (bsdv[i] == 1)
               {
                  Ai = A + ldA*i;
                  Work[i] = Ai[i-1];
                  Ai[i-1] = zero;
                  --i;
               }
               else
                  Work[i] = zero;
               --i;
            }
            
            i = n-1;
            while (i >= 0)
            {
               if (bsdv[i] == 0)
               {
                  ip = ipiv[i];
                  if (ip != i && i < n-1)
                  {
                     for (int_t j = i+1; j < n; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i];
                        Aj[i] = temp;
                     }
                  }
               }
               else
               {
                  ip = ipiv[i];
                  if (ip != i && i < n-1)
                  {
                     for (int_t j = i+1; j < n; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i-1];
                        Aj[i-1] = temp;
                     }
                  }
                  --i;
               }
               --i;
            }
         }
         else
         {
            i = 0;
            while (i < n)
            {
               if (bsdv[i] == 0)
               {
                  ip = ipiv[i];
                  if (ip != i && i < n-1)
                  {
                     for (int_t j = i+1; j < n; ++j)
                     {
                        real_t * Aj = A+ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i];
                        Aj[i] = temp;
                     }
                  }
               }
               else
               {
                  ip = ipiv[i];
                  ++i;
                  if (ip != i && i < n-1)
                  {
                     for (int_t j = i+1; j < n; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i-1];
                        Aj[i-1] = temp;
                     }
                  }
               }
               ++i;
            }
            i = n-1;
            while (i >= 0)
            {
               if (bsdv[i] == 1)
               {
                  real_t * Ai = A + ldA*i;
                  Ai[i-1] = Work[i];
                  --i;
               }
               --i;
            }
         }
      }
      else //lower
      {
         if (way == 'C' || way == 'c')
         {
            real_t * Ai;
            Work[n-1] = zero;
            while (i < n)
            {
               if ( bsdv[i] == 1 && i < n-1)
               {
                  Ai = A + ldA*i;
                  Work[i] = Ai[i+1];
                  Ai[i+1] = zero;
                  ++i;
               }
               Work[i] = zero;
               ++i;
            }
            i = 0;
            while (i < n)
            {
               if (bsdv[i] == 0)
               {
                  ip = ipiv[i];
                  if (ip != i && i > 0)
                  {
                     for (int_t j = 0; j < i; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i];
                        Aj[i] = temp;
                     }
                  }
               }
               else
               {
                  ip = ipiv[i];
                  if (ip != i+1 && i > 0)
                  {
                     for (int_t j = 0; j < i; ++j)
                     {
                        real_t * Aj = A + ldA*j;
                        real_t temp = Aj[ip];
                        Aj[ip] = Aj[i+1];
                        Aj[i+1] = temp;
                     }
                  }
                  ++i;
               }
               ++i;
            }
         }
         else
         {
            i = n-1;
            while (i >= 0)
            {
               if (bsdv[i] == 0)
               {
                  ip = ipiv[i];
                  if (ip != i && i > 0)
                  {
                     for (int_t j = 0; j < i; ++j)
                     {
                        real_t * Aj = A+ldA*j;
                        real_t temp = Aj[i];
                        Aj[i] = Aj[ip];
                        Aj[ip] = temp;
                     }
                  }
               }
               else
               {
                  ip = ipiv[i];
                  --i;
                  if (ip != i && i > 0)
                  {
                     for (int_t j = 0; j < i; ++j)
                     {
                        real_t * Aj = A+ldA*j;
                        real_t temp = Aj[i+1];
                        Aj[i+1] = Aj[ip];
                        Aj[ip] = temp;
                     }
                  }
               }
               --i;
            }
            i = 0;
            while (i < n-1)
            {
               if (bsdv[i] == 1)
               {
                  real_t * Ai = A+ldA*i;
                  Ai[i+1] = Work[i];
                  ++i;
               }
               ++i;
            }
         }
      }
      return 0;
   }
   
   /// @brief Converts A, given by sytrf or hetrf, into L and D and vice-versa.
   ///
   /// The off-diagonal elements of D are returned in Work.  Also applies (or reverses) permutation done in TRF.
   ///
   /// This is an auxiliary routine for sytrs and hetrs.
   /// @return 0
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric or Hermitian matrix A is stored as upper triangular or lower triangular.
   /// @param way Indicates whether the matrix A should be converted to L and D, or if the changes need to be reverted.
   ///
   ///      'C' = convert
   ///      'R' = revert
   ///
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Complex symmetric matrix size ldA-by-n.  On entry, the block diagonal matrix D and the multipliers used to obtain the factor U or L as computed by sytrf or hetrf.
   /// @param ldA Column length or the matrix A.  ldA >= n
   /// @param ipiv Integer array size n.  On entry, should contain the details of the interchanges of D.
   /// @param bsdv Bool array size n.  On entry, should contain the details of the block structure of D, as computed by sytrf or hetrf.
   /// @param Work Real array length n.  On exit, contains the off-diagonal elements of A, with any other elements set to zero.
   /// @ingroup AUX

   template< typename real_t>
   int_t syconv(const char uplo, const char way, const int_t n, complex<real_t> * const A, const int_t ldA, int_t * ipiv, bool * bsdv, real_t * Work)
   {
      return latl::syconv(uplo, way, n, A, ldA, ipiv, bsdv, Work);
   }
}


#endif
