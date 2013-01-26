//
//  lange.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 6/1/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lange_h
#define _lange_h

/// @file lange.h Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a general matrix A.

#include "lassq.h"
#include <cmath>
#include "latl.h"

namespace latl
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value
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
   /// @param m Number of rows to be included in the norm. m >= 0
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param A Real matrix size m-by-n.
   /// @param ldA Column length of the matrix A.  ldA >= m
   /// @ingroup NORM

   template< typename real_t>
   real_t lange(const char normType, const int_t m, const int_t n, real_t * const A, const int_t ldA)
   {
      using std::isnan;
      using std::abs;
      using std::sqrt;
      const real_t zero(0.0);
      real_t value(0.0);
      if (m == 0 || n == 0)
         return zero;
      if (normType == 'M' || normType == 'm')
      {
         real_t * Aj = A;
         real_t temp;
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = 0; i < m; ++i)
            {
               temp = abs(Aj[i]);
               if (temp > value)
               {
                  value = temp;
               }
               else if (isnan(temp))
                  return temp;
            }
            Aj += ldA;
         }
      }
      else if ( normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t sum, * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            sum = zero;
            for (int_t i = 0; i < m; ++i)
            {
               sum += abs(Aj[i]);
            }
            if (sum > value)
            {
               value = sum;
            }
            else if (isnan(sum))
               return sum;
            Aj += ldA;
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         real_t * Work = new real_t[m];
         real_t temp;
         for (int_t i = 0; i < m; ++i)
         {
            Work[i] = abs(A[i]);
         }
         real_t * Aj = A+ldA;
         for (int_t j = 1; j < n; ++j)
         {
            for (int_t i = 0; i < m; ++i)
            {
               Work[i] += abs(Aj[i]);
            }
            Aj += ldA;
         }
         for (int_t i = 0; i < m; ++i)
         {
            temp = Work[i];
            if (temp > value)
            {
               value = temp;
            }
            else if (isnan(temp))
            {
               delete [] Work;
               return temp;
            }
         }
         delete [] Work;
      }
      else if ( normType == 'F' || normType == 'f' || normType == 'E' ||normType == 'e')
      {
         real_t scale(0.0), sum(1.0), * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            latl::lassq(m, Aj, 1, scale, sum);
            Aj += ldA;
         }
         value = scale*sqrt(sum);
      }
      return value;
   }
   
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value
   ///
   /// @tparam real_t Real floating point type
   /// @return Calculated norm value for the specified type.
   /// @param normType Type should be specified as follows:
   ///
   ///     'M' or 'm' = maximum absolute value over all elements in A.
   ///         Note: this is not a consistent matrix norm.
   ///     '1', 'O', or 'o' = one norm of the matrix A, the maximum value of the sums of each column.
   ///     'I' or 'i' = the infinity norm of the matrix A, the maximum value of the sum of each row.
   ///     'F', 'f', 'E', or 'e' = the Frobenius norm of the matrix A.
   ///         This the square root of the sum of the squares of each element in A.
   ///
   /// @param m Number of rows to be included in the norm. m >= 0
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param A Real matrix size m-by-n.
   /// @param ldA Column length of the matrix A.  ldA >= m
   /// @ingroup NORM
   
   template< typename real_t>
   real_t lange(const char normType, const int_t m, const int_t n, complex<real_t> * const A, const int_t ldA)
   {
      using std::isnan;
      using std::abs;
      const real_t zero(0.0);
      real_t value(0.0);
      if (m == 0 || n == 0)
         return zero;
      if (normType == 'M' || normType == 'm')
      {
         complex< real_t> * Aj = A;
         real_t temp;
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = 0; i < m; ++i)
            {
               temp = abs(Aj[i]);
               if (temp > value)
               {
                  value = temp;
               }
               else if (isnan(temp))
                  return temp;
            }
            Aj += ldA;
         }
      }
      else if ( normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t sum;
         complex<real_t> * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            sum = zero;
            for (int_t i = 0; i < m; ++i)
            {
               sum += abs(Aj[i]);
            }
            if (sum > value)
            {
               value = sum;
            }
            else if (isnan(sum))
               return sum;
            Aj += ldA;
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         real_t * Work = new real_t[m];
         real_t temp;
         for (int_t i = 0; i < m; ++i)
         {
            Work[i] = abs(A[i]);
         }
         complex< real_t> * Aj = A+ldA;
         for (int_t j = 1; j < n; ++j)
         {
            for (int_t i = 0; i < m; ++i)
            {
               Work[i] += abs(Aj[i]);
            }
            Aj += ldA;
         }
         for (int_t i = 0; i < m; ++i)
         {
            temp = Work[i];
            if (temp > value)
            {
               value = temp;
            }
            else if (isnan(temp))
            {
               delete [] Work;
               return temp;
            }
         }
         delete [] Work;
      }
      else if ( normType == 'F' || normType == 'f' || normType == 'E' ||normType == 'e')
      {
         real_t scale(0.0), sum(1.0);
         complex< real_t> * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            latl::lassq(m, Aj, 1, scale, sum);
            Aj += ldA;
         }
         value = scale*sqrt(sum);
      }
      return value;
   }
}
#endif
