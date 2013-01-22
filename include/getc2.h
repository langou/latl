//
//  getc2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 6/21/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _getc2_h
#define _getc2_h

/// @file getc2.h Computes an LU factorization with complete pivoting of an n-by-n matrix.

#include <limits>
#include "labad.h"
#include "swap.h"
#include "ger.h"
#include "latl.h"

namespace latl
{
   /// @brief Computes an LU factorization with complete pivoting of an n-by-n matrix.
   ///
   /// The factorization has the form
   /// A = P * L * U * Q
   /// where P and Q are permutation matrices.  L is lower triangular with unit diagonal elements and U is upper triangular.
   ///
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Real array size ldA-by-n.  On entry, the n-by-n matrix to be factored.  On exit, contains the factors L and U from the factorization.  The unit diagonal elements are not stored.  If U(k,k) appears to be less than smin, U(k,k) is given the value of smin, giving a nonsingular perturbed system.
   /// @param ldA The column length of the array A.  ldA >= n
   /// @param IPIV Integer array size n.  Contains pivot indices; for 0 <= i <= n-1, row i of the matrix has been exchanged with row IPIV[i].
   /// @param JPIV Integer array size n.  Contains pivot indices; for 0 <= j <= n-1, column j of the matrix has been exchanged with column JPIV[j].
   /// @ingroup TRF
   
   template<typename real_t>
   int getc2(const int_t n, real_t * const A, const int_t ldA, int_t * const IPIV, int_t * const JPIV)
   {
      using std::numeric_limits;

      if (n < 0)
         return -1;
      if (ldA < n)
         return -3;
      
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t eps = numeric_limits<real_t>::epsilon();
      real_t smlnum = numeric_limits<real_t>::min()/eps;
      real_t bignum = one/smlnum;
      latl::labad(smlnum, bignum);
      
      int_t ipv, jpv, info = 0;
      real_t xmax, smin;
      real_t * Ai = A;
      real_t * Aij = A;
      for (int_t i = 0; i < n-1; ++i)
      {
         xmax = zero;
         for (int_t ip = i; ip < n; ++ip)
         {
            Ai = A + ip;
            for (int_t jp = i; jp < n; ++jp)
            {
               Aij = Ai + ldA*jp;
               if (abs(Aij[0]) >= xmax)
               {
                  xmax = abs(Aij[0]);
                  ipv = ip;
                  jpv = jp;
               }
            }
         }
         if ( i == 0)
         {
            smin = std::max(eps*xmax, smlnum);
         }
         
         if (ipv != i)
         {
            latl::swap(n, A+ipv, ldA, A+i, ldA);
         }
         IPIV[i] = ipv;
         
         if (jpv != i)
         {
            latl::swap(n, A+ldA*jpv, 1, A+ldA*i, 1);
         }
         JPIV[i] = jpv;
         
         Aij = A+i+ldA*i;
         if (abs(Aij[0]) < smin)
         {
            if (info == 0)
               info = i+1;
            Aij[0] = smin;
         }
         
         for (int_t k = 1; k < (n-i); ++k)
         {
            Aij[k] = Aij[k] / Aij[0];
         }
         latl::ger(n-i-1, n-i-1, -one, Aij+1, 1, Aij+ldA, ldA, Aij+1+ldA, ldA);
      }
      
      Aij = A+ldA*(n-1)+(n-1);
      if (abs(Aij[0]) < smin)
      {
         if (info == 0)
            info = n;
         Aij[0] = smin;
      }
      JPIV[n-1] = IPIV[n-1] = n-1;
      
      return info;
   }
   
   /// @brief Computes an LU factorization with complete pivoting of an n-by-n matrix.
   ///
   /// The factorization has the form
   /// A = P * L * U * Q
   /// where P and Q are permutation matrices.  L is lower triangular with unit diagonal elements and U is upper triangular.
   ///
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Complex array size ldA-by-n.  On entry, the n-by-n matrix to be factored.  On exit, contains the factors L and U from the factorization.  The unit diagonal elements are not stored.  If U(k,k) appears to be less than smin, U(k,k) is given the value of smin, giving a nonsingular perturbed system.
   /// @param ldA The column length of the array A.  ldA >= n
   /// @param IPIV Integer array size n.  Contains pivot indices; for 0 <= i <= n-1, row i of the matrix has been exchanged with row IPIV[i].
   /// @param JPIV Integer array size n.  Contains pivot indices; for 0 <= j <= n-1, column j of the matrix has been exchanged with column JPIV[j].
   /// @ingroup TRF
   
   template<typename real_t>
   int getc2(const int_t n, complex<real_t> * const A, const int_t ldA, int_t * const IPIV, int_t * const JPIV)
   {
      using std::numeric_limits;

      if (n < 0)
         return -1;
      if (ldA < n)
         return -3;
      
      const real_t one(1.0);
      const complex<real_t> cxone(1.0);
      const real_t zero(0.0);
      const real_t eps = numeric_limits<real_t>::epsilon();
      real_t smlnum = numeric_limits<real_t>::min()/eps;
      real_t bignum = one/smlnum;
      latl::labad(smlnum, bignum);
      
      int_t ipv, jpv, info = 0;
      real_t xmax;
      complex<real_t> cmin;
      complex<real_t> * Ai = A;
      complex<real_t> * Aij = A;
      for (int_t i = 0; i < n-1; ++i)
      {
         xmax = zero;
         for (int_t ip = i; ip < n; ++ip)
         {
            Ai = A + ip;
            for (int_t jp = i; jp < n; ++jp)
            {
               Aij = Ai + ldA*jp;
               if (abs(Aij[0]) >= xmax)
               {
                  xmax = abs(Aij[0]);
                  ipv = ip;
                  jpv = jp;
               }
            }
         }
         if ( i == 0)
         {
            cmin = std::max(eps*xmax, smlnum);
         }
         
         if (ipv != i)
         {
            latl::swap(n, A+ipv, ldA, A+i, ldA);
         }
         IPIV[i] = ipv;
         
         if (jpv != i)
         {
            latl::swap(n, A+ldA*jpv, 1, A+ldA*i, 1);
         }
         JPIV[i] = jpv;
         
         Aij = A+i+ldA*i;
         if (abs(Aij[0]) < real(cmin))
         {
            if (info == 0)
               info = i+1;
            Aij[0] = cmin;
         }
         
         for (int_t k = 1; k < (n-i); ++k)
         {
            Aij[k] = Aij[k] / Aij[0];
         }
         latl::ger(n-i-1, n-i-1, -cxone, Aij+1, 1, Aij+ldA, ldA, Aij+1+ldA, ldA);
      }
      
      ipv = (n-1)*(n-1);
      if (abs(A[ipv]) < real(cmin))
      {
         if (info == 0)
            info = n;
         A[ipv] = cmin;
      }
      JPIV[n-1] = IPIV[n-1] = n-1;
      
      return info;
   }
}

#endif
