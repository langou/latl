//
//  geequ.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 6/2/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _geequ_h
#define _geequ_h

/// @file geequ.h Computes row and column scalings intended to equilibrate and reduce the condition number of a general m-by-n matrix.

#include <limits>
#include <cmath>
#include <algorithm>
#include "latl.h"

namespace LATL
{
   /// @brief Computes row and column scalings intended to equilibrate an m-by-n matrix A and reduce its condition number.
   ///
   /// The row scale factors R and and column scale factors C are chosen to try to make the largest element in each row and column of the matrix B, with B(i, j) = R[i] * A(i,j) * C[j], have absolute value 1.
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith row of A is exactly 0.
   /// @return m+j+1 if the jth column of A is exactly 0.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of the columns of the matrix A. n >= 0
   /// @param A Real matrix size m-by-n whose equilibration factors are to be computed.
   /// @param ldA Column length of A.  ldA >= m
   /// @param R Real array of size m.  On exit, if the return value is 0 or is greater than m, R contains the row scale factors for A.  Scale factors are restricted to values between the smallest safe number (smlnum) and largest safe number (bignum).
   /// @param C Real array of size n.  On exit, if the return value equals 0, C contains the column scale factors for A.  Scale factors are restricted to values between smlnum and bignum.
   /// @param rowcnd Real number.  On exit, if the return value is 0 or is greater than m, rowcnd is the ratio of the smallest R(i) to the largest R(i).  If rowcnd >= .1 and amax is neither too large nor too small, it is not worth scaling by R.
   /// @param colcnd Real number.  On exit, if the return value is 0, colcnd is the ratio of the smallest C(i) to the largest C(i).  If colcnd >= .1, it is not worth scaling by C.
   /// @param amax Real number.  On exit, the absolute value of the largest matrix element.  If amax if very close to overflow or underflow, the matrix should be scaled.
   /// @ingroup COMP
   
   template< typename real_t>
   int GEEEQU(const int_t m, const int_t n, real_t * const A, const int_t ldA, real_t * const R, real_t * const C, real_t &rowcnd, real_t &colcnd, real_t &amax)
   {
      using std::numeric_limits;
      
      if ( m < 0)
         return -1;
      if ( n < 0)
         return -2;
      if (ldA < m)
         return -4;
      
      if ( m == 0 || n == 0)
      {
         rowcnd = 1;
         colcnd = 1;
         amax = 0;
         return 0;
      }
      
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t smlnum = numeric_limits<real_t>::min();
      const real_t bignum = one/smlnum;
      
      for (int_t i = 0; i < m; ++i)
      {
         R[i] = std::abs(A[i]);
      }
      real_t * Aj = A + ldA;
      for (int_t j = 1; j < n; ++j)
      {
         for (int_t i = 0; i < m; ++i)
         {
            if (std::abs(Aj[i]) > R[i])
            {
               R[i] = std::abs(Aj[i]);
            }
         }
         Aj += ldA;
      }
      real_t rcmin = bignum;
      real_t rcmax = zero;
      for (int_t i = 0; i < m; ++i)
      {
         if (R[i] > rcmax)
            rcmax = R[i];
         if (R[i] < rcmin)
            rcmin = R[i];
      }
      amax = rcmax;
      
      if (rcmin == zero)
      {
         for (int_t i = 0; i < m; ++i)
         {
            if (R[i] == zero)
            {
               return i+1;
            }
         }
      }
      
      for (int_t i = 0; i < m; ++i)
      {
         R[i] = one/std::min(std::max(R[i], smlnum), bignum);
      }
      
      rowcnd = std::max(rcmin, smlnum)/ std::min(rcmax, bignum);
      
      Aj = A;
      for (int_t j = 0; j < n; ++j)
      {
         C[j] = std::abs(Aj[0]*R[0]);
         Aj += ldA;
      }
      Aj = A;
      for (int_t j = 0; j < n; ++j)
      {
         for (int_t i = 1; i < m; ++i)
         {
            C[j] = std::max(C[j], std::abs(Aj[i]*R[i]));
         }
         Aj += ldA;
      }
      
      rcmin = bignum;
      rcmax = zero;
      for (int_t j = 0; j < n; ++j)
      {
         if (C[j] < rcmin)
         {
            rcmin = C[j];
         }
         if (C[j] > rcmax)
         {
            rcmax = C[j];
         }
      }
      if (rcmin == zero)
      {
         for (int_t j = 0; j < n; ++j)
         {
            if (C[j] == zero)
            {
               return m+j+1;
            }
         }
      }
      
      for (int_t j = 0; j < n; ++j)
      {
         C[j] = one/std::min(std::max(C[j], smlnum), bignum);
      }
      
      colcnd = std::max(rcmin, smlnum)/std::min(rcmax, bignum);
      
      return 0;
   }
   
   /// @brief Computes row and column scalings intended to equilibrate an m-by-n matrix A and reduce its condition number.
   ///
   /// The row and column scale factors are chosen to try to make the largest element in each row and column of the matrix B, with B(i, j) = R[i] * A(i,j) * C[j], have absolute value 1.
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith row of A is exactly 0.
   /// @return m+j+1 if the jth column of A is exactly 0.
   /// @tparam real_t Floating point type.
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of the columns of the matrix A. n >= 0
   /// @param A Complex matrix size m-by-n whose equilibration factors are to be computed.
   /// @param ldA Column length of A.  ldA >= m
   /// @param R Real array of size m.  On exit, if the return value is 0 or is greater than m, R contains the row scale factors for A.  Restricted to values between the smallest safe number (smlnum) and largest safe number (bignum).
   /// @param C Real array of size n.  On exit, if the return value equals 0, C contains the column scale factors for A.  Restricted to values between smlnum and bignum.
   /// @param rowcnd Real number.  On exit, if the return value is 0 or is greater than m, rowcnd is the ratio of the smallest R(i) to the largest R(i).  If rowcnd >= .1 and amax is neither too large nor too small, it is not worth scaling by R.
   /// @param colcnd Real number.  On exit, if the return value is 0, colcnd is the ratio of the smallest C(i) to the largest C(i).  If colcnd >= .1, it is not worth scaling by C.
   /// @param amax Real number.  On exit, the absolute value of the largest matrix element.  If amax if very close to overflow or underflow, the matrix should be scaled.
   /// @ingroup COMP

   template< typename real_t>
   int GEEEQU(const int_t m, const int_t n, complex<real_t> * const A, const int_t ldA, real_t * const R, real_t * const C, real_t &rowcnd, real_t &colcnd, real_t &amax)
   {
      using std::numeric_limits;

      if ( m < 0)
         return -1;
      if ( n < 0)
         return -2;
      if (ldA < m)
         return -4;
      
      if ( m == 0 || n == 0)
      {
         rowcnd = 1;
         colcnd = 1;
         amax = 0;
         return 0;
      }
      
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t smlnum = numeric_limits<real_t>::min();
      const real_t bignum = one/smlnum;
      
      for (int_t i = 0; i < m; ++i)
      {
         R[i] = std::abs(real(A[i]))+ std::abs(imag(A[i]));
      }
      complex<real_t> * Aj = A + ldA;
      for (int_t j = 1; j < n; ++j)
      {
         for (int_t i = 0; i < m; ++i)
         {
            R[i] = std::max(R[i], std::abs(real(Aj[i]))+ std::abs(imag(Aj[i])));
         }
         Aj += ldA;
      }
      real_t rcmin = bignum;
      real_t rcmax = zero;
      for (int_t i = 0; i < m; ++i)
      {
         if (R[i] > rcmax)
            rcmax = R[i];
         if (R[i] < rcmin)
            rcmin = R[i];
      }
      amax = rcmax;
      
      if (rcmin == zero)
      {
         for (int_t i = 0; i < m; ++i)
         {
            if (R[i] == zero)
            {
               return i+1;
            }
         }
      }
      
      for (int_t i = 0; i < m; ++i)
      {
         R[i] = one/std::min(std::max(R[i], smlnum), bignum);
      }
      
      rowcnd = std::max(rcmin, smlnum)/ std::min(rcmax, bignum);
      
      Aj = A;
      for (int_t j = 0; j < n; ++j)
      {
         C[j] = (std::abs(real(Aj[0])) + std::abs(imag(Aj[0])))*R[0];
         Aj += ldA;
      }
      Aj = A;
      for (int_t j = 0; j < n; ++j)
      {
         for (int_t i = 1; i < m; ++i)
         {
            C[j] = std::max(C[j], (std::abs(real(Aj[i])) + std::abs(imag(Aj[i])))*R[i]);
         }
         Aj += ldA;
      }
      
      rcmin = bignum;
      rcmax = zero;
      for (int_t j = 0; j < n; ++j)
      {
         if (C[j] < rcmin)
         {
            rcmin = C[j];
         }
         if (C[j] > rcmax)
         {
            rcmax = C[j];
         }
      }
      if (rcmin == zero)
      {
         for (int_t j = 0; j < n; ++j)
         {
            if (C[j] == zero)
            {
               return m+j+1;
            }
         }
      }
      
      for (int_t j = 0; j < n; ++j)
      {
         C[j] = one/std::min(std::max(C[j], smlnum), bignum);
      }
      
      colcnd = std::max(rcmin, smlnum)/std::min(rcmax, bignum);
      
      return 0;
      
   }
   
}

#endif
