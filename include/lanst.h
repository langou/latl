//
//  lanst.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/20/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lanst_h
#define _lanst_h

/// @file lanst.h Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a symmetric tridiagonal matrix A.

#include "latl.h"
#include <cmath>
#include "lassq.h"

namespace latl
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a real symmetric tridiagonal matrix.
   ///
   /// @return Calculated norm value for the specified type.
   /// @tparam real_t Real floating point type
   /// @param normType Type should be specified as follows:
   ///
   ///     'M' or 'm' = maximum absolute value over all elements in A.
   ///         Note: this is not a consistent matrix norm.
   ///     '1', 'O', or 'o' = one norm of the matrix A, the maximum value of the sums of each column.
   ///     'I' or 'i' = the infinity norm of the matrix A, the maximum value of the sum of each row.
   ///     'F', 'f', 'E', or 'e' = the Frobenius norm of the matrix A.
   ///         This the square root of the sum of the squares of each element in A.
   ///
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param D Real array, size n.  On entry, D should contain the n diagonal elements of A.
   /// @param E Real array, size (n-1).  On entry, E should contain the (n-1) superdiagonal (or subdiagonal) elements of A.
   /// @ingroup NORM
   
   template <typename real_t>
   real_t lanst(const char normType, const int_t n, real_t * const D, real_t * const E)
   {
      using std::isnan;
      using std::abs;
      real_t value(0.0);
      if (n <= 0)
         return value;
      if (normType == 'M' || normType == 'm')
      {
         value = abs(D[n-1]);
         real_t temp;
         for (int_t i = 0; i < n-1; ++i)
         {
            temp = abs(D[i]);
            if (temp > value)
               value = temp;
            else if (isnan(temp))
               return temp;
            temp = abs(E[i]);
            if (temp > value)
               value = temp;
            else if (isnan(temp))
               return temp;
         }
      }
      else if (normType == 'O' || normType == 'o' || normType == '1' || normType == 'I' || normType == 'i')
      {
         if (n == 1)
            value = abs(D[0]);
         else
         {
            real_t sum = abs(E[n-2]) + abs(D[n-1]);
            value = abs(D[0]) + abs(E[0]);
            if (value < sum)
               value = sum;
            else if (isnan(sum))
               return sum;
            for (int_t i = 1; i < n-1; ++i)
            {
               sum = abs(D[i]) + abs(E[i]) + abs(E[i-1]);
               if (value < sum)
               {
                  value = sum;
               }
               else if (isnan(sum))
                  return sum;
            }
         }
      }
      else if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale(0.0), sum(1.0);
         if (n > 1)
         {
            latl::lassq(n-1, E, 1, scale, sum);
            sum *= 2;
         }
         latl::lassq(n, D, 1, scale, sum);
         value = scale*sqrt(sum);
      }
      return value;
   }
   
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a complex symmetric tridiagonal matrix.
   ///
   /// @return Calculated norm value for the specified type.
   /// @tparam real_t Real floating point type
   /// @param normType Type should be specified as follows:
   ///
   ///     'M' or 'm' = maximum absolute value over all elements in A.
   ///         Note: this is not a consistent matrix norm.
   ///     '1', 'O', or 'o' = one norm of the matrix A, the maximum value of the sums of each column.
   ///     'I' or 'i' = the infinity norm of the matrix A, the maximum value of the sum of each row.
   ///     'F', 'f', 'E', or 'e' = the Frobenius norm of the matrix A.
   ///         This the square root of the sum of the squares of each element in A.
   ///
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param D Complex array, size n.  On entry, D should contain the n diagonal elements of A.
   /// @param E Complex array, size (n-1).  On entry, E should contain the (n-1) superdiagonal (or subdiagonal) elements of A.
   /// @ingroup NORM
   
   template <typename real_t>
   real_t lanst(const char normType, const int_t n, complex<real_t> * const D, complex<real_t> * const E)
   {
      using std::isnan;
      using std::abs;
      real_t value(0.0);
      if (n <= 0)
         return value;
      if (normType == 'M' || normType == 'm')
      {
         value = abs(D[n-1]);
         real_t temp;
         for (int_t i = 0; i < n-1; ++i)
         {
            temp = abs(D[i]);
            if (temp > value)
               value = temp;
            else if (isnan(temp))
               return temp;
            temp = abs(E[i]);
            if (temp > value)
               value = temp;
            else if (isnan(temp))
               return temp;
         }
      }
      else if (normType == 'O' || normType == 'o' || normType == '1' || normType == 'I' || normType == 'i')
      {
         if (n == 1)
            value = abs(D[0]);
         else
         {
            real_t sum = abs(E[n-2]) + abs(D[n-1]);
            value = abs(D[0]) + abs(E[0]);
            if (value < sum)
               value = sum;
            else if (isnan(sum))
               return sum;
            for (int_t i = 1; i < n-1; ++i)
            {
               sum = abs(D[i]) + abs(E[i]) + abs(E[i-1]);
               if (value < sum)
               {
                  value = sum;
               }
               else if (isnan(sum))
                  return sum;
            }
         }
      }
      else if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale(0.0), sum(1.0);
         if (n > 1)
         {
            latl::lassq(n-1, E, 1, scale, sum);
            sum *= 2;
         }
         latl::lassq(n, D, 1, scale, sum);
         value = scale*sqrt(sum);
      }
      return value;
   }
}
#endif
