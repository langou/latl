//
//  getf2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 5/17/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _getf2_h
#define _getf2_h

/// @file getf2.h Computes an LU factorization of an m-by-n matrix A using partial pivoting.

#include <algorithm>
#include <cmath>
#include <limits>
#include "imax.h"
#include "swap.h"
#include "scal.h"
#include "ger.h"
#include "latl.h"

namespace latl
{
   /// @brief Computes an LU factorization of an m-by-n matrix A using partial pivoting.
   ///
   /// The factorization has the form
   ///
   ///     A = P * L * U
   ///
   /// where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower trapezoidal if m > n), and U is upper triangular (upper trapezoidal if m < n).
   /// This is the right-looking Level 2 BLAS version of the algorithm.
   /// @return  0 if success
   /// @return -i if the ith argument is invalid.
   /// @return  i+1 if the ith column of matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param A Real matrix size m-by-n.  On exit, the factors L and U.
   /// @param ldA Column length of matrix A.  ldA >= m
   /// @param IPIV Permutation matrix size min(m,n).  On exit, row k of A was exchanged with IPIV[k].
   /// @ingroup TRF
   
   template< typename real_t>
   int_t getf2(const int_t m, const int_t n, real_t * const A, const int_t ldA,  int_t *const IPIV)
   {
      using std::numeric_limits;

      if (m < 0)
         return -1;
      if (n < 0)
         return -2;
      if (ldA < m)
         return -4;
      
      if (m == 0 || n == 0)
         return 0;
      
      const real_t sfmin = numeric_limits<real_t>::min();
      const real_t one = 1.0;
      const real_t zero = 0.0;
      
      int_t jp = 0, info = 0;
      real_t *Ajj = A;
      real_t *Ajjp1 = Ajj+1;
      
      for (int_t j = 0; j < std::min(m,n); ++j)
      {
         jp = j + latl::imax(m-j, Ajj, 1);
         IPIV[j]=jp;
         if (A[ldA*j+jp] != zero)
         {
            if (jp != j)
            {
               latl::swap(n, A+j, ldA, A+jp, ldA);
            }
            if (j < m)
            {
               if (std::abs(Ajj[0]) >= sfmin)
               {
                  latl::scal(m-j-1, one/Ajj[0], Ajjp1, 1);
               }
               else
               {
                  for (int_t i = 1; i < m-j; ++i)
                  {
                     Ajj[i] = Ajj[i]/Ajj[0];
                  }
               }
            }
         }
         else
         {
            if (info == 0)
            {
               info = j+1;
            }
         }
         if (j < std::min(m,n))
         {
            latl::ger(m-j-1, n-j-1, -one, Ajjp1, 1, Ajj+ldA, ldA, Ajjp1+ldA, ldA);
         }
         Ajj+= (ldA+1);
         Ajjp1 = Ajj+1;
      }
      
      return info;
   }
   
   /// @brief Computes an LU factorization of an m-by-n matrix A using partial pivoting
   ///
   /// The factorization has the form
   ///
   ///     A = P * L * U
   ///
   /// where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower trapezoidal if m > n), and U is upper triangular (upper trapezoidal if m < n).
   /// This is the right-looking Level 2 BLAS version of the algorithm.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column of matrix A has a zero pivot, meaning U is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param A Complex matrix size m-by-n.  On exit, the factors L and U.
   /// @param ldA Column length of matrix A.  ldA >= m
   /// @param IPIV Permutation matrix size min(m,n).  On exit, row k of A was exchanged with IPIV[k].
   /// @ingroup TRF
   
   template< typename real_t>
   int_t getf2(const int_t m, const int_t n, complex<real_t> * const A, const int_t ldA,  int_t * const IPIV)
   {
      using std::numeric_limits;

      if (m < 0)
         return -1;
      if (n < 0)
         return -2;
      if (ldA < m)
         return -4;
      
      if (m == 0 || n == 0)
         return 0;
      
      const real_t sfmin = numeric_limits<real_t>::min();
      const complex<real_t> one(1.0, 0);
      const complex<real_t> zero(0.0,0.0);
      
      int_t jp = 0, info = 0;
      complex<real_t> *Ajj = A;
      complex<real_t> *Ajjp1 = Ajj+1;
      
      for (int_t j = 0; j < std::min(m,n); ++j)
      {
         jp = j + latl::imax(m-j, Ajj, 1);
         IPIV[j]=jp;
         if (A[ldA*j+jp] != zero)
         {
            if (jp != j)
            {
               latl::swap(n, A+j, ldA, A+jp, ldA);
            }
            if (j < m)
            {
               if (std::abs(Ajj[0]) >= sfmin)
               {
                  latl::scal(m-j-1, one/Ajj[0], Ajjp1, 1);
               }
               else
               {
                  for (int_t i = 1; i < m-j; ++i)
                  {
                     Ajj[i] = Ajj[i]/Ajj[0];
                  }
               }
            }
         }
         else
         {
            if (info == 0)
            {
               info = j+1;
            }
         }
         if (j < std::min(m,n))
         {
            latl::ger(m-j-1, n-j-1, -one, Ajjp1, 1, Ajj+ldA, ldA, Ajjp1+ldA, ldA);
         }
         Ajj+= (ldA+1);
         Ajjp1 = Ajj+1;
      }
      return info;
   }
   
}
#endif
