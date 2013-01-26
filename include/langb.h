//
//  langb.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/19/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _langb_h
#define _langb_h

/// @file langb.h Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a banded matrix A.

#include "latl.h"
#include <cmath>
#include "lassq.h"

namespace latl
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a real banded matrix.
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
   /// @param kL The number of subdiagonals within the band of A.  kL >= 0.
   /// @param kU The number of superdiagonals within the band of A.  kU >= 0.
   /// @param AB Real matrix size ldAB-by-n.  On entry, the matrix A in band storage.
   /// @param ldAB Column length of the matrix AB.  ldAB >= kU+kL+1
   /// @ingroup NORM
   
   template< typename real_t>
   real_t langb(const char normType, const int_t n, const int_t kL, const int_t kU, real_t * const AB, const int_t ldAB)
   {
      using std::isnan;
      using std::abs;
      using std::min;
      using std::max;
      const real_t zero(0.0);
      const int_t intzero(0);
      real_t value(0.0);
      if (n == 0)
         return value;
      if (normType == 'M' || normType == 'm')
      {
         real_t * ABj = AB;
         real_t temp;
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = max(kU-j, intzero); i < min(n+kU-j, kL+kU+1); ++i)
            {
               temp = abs(ABj[i]);
               if (temp > value)
               {
                  value = temp;
               }
               else if (isnan(temp))
                  return temp;
            }
            ABj += ldAB;
         }
      }
      else if ( normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t sum, * ABj = AB;
         for (int_t j = 0; j < n; ++j)
         {
            sum = zero;
            for (int_t i = max(kU-j, intzero); i < min(n+kU-j, kL+kU+1); ++i)
            {
               sum += abs(ABj[i]);
            }
            if (sum > value)
            {
               value = sum;
            }
            else if (isnan(sum))
               return sum;
            ABj += ldAB;
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         real_t * Work = new real_t[n];
         real_t temp;
         for (int_t i = 0; i < n; ++i)
         {
            Work[i] = zero;
         }
         real_t * ABj = AB;
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = max(kU-j, intzero); i < min(n+kU-j, kL+kU+1); ++i)
            {
               Work[j+i-kU] += abs(ABj[i]);
            }
            ABj += ldAB;
         }
         for (int_t i = 0; i < n; ++i)
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
         real_t scale(0.0), sum(1.0), * ABj = AB;
         int_t L, K;
         for (int_t j = 0; j < n; ++j)
         {
            L = max(intzero,j-kU);
            K = kU -j+L;
            latl::lassq(min(n, j+kL+1)-L, ABj+K, 1, scale, sum);
            ABj += ldAB;
         }
         value = scale*sqrt(sum);
      }
      return value;
   }

   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a complex banded matrix.
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
   /// @param kL The number of subdiagonals within the band of A.  kL >= 0.
   /// @param kU The number of superdiagonals within the band of A.  kU >= 0.
   /// @param AB Real matrix size ldAB-by-n.  On entry, the matrix A in band storage.
   /// @param ldAB Column length of the matrix AB.  ldAB >= kL+kU+1
   /// @ingroup NORM
   
   template< typename real_t>
   real_t langb(const char normType, const int_t n, const int_t kL, const int_t kU, complex<real_t> * const AB, const int_t ldAB)
   {
      using std::isnan;
      using std::abs;
      using std::min;
      using std::max;
      const real_t zero(0.0);
      const int_t intzero(0);
      real_t value(0.0);
      if (n == 0)
         return value;
      if (normType == 'M' || normType == 'm')
      {
         complex<real_t> * ABj = AB;
         real_t temp;
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = max(kU-j, intzero); i < min(n+kU-j, kL+kU+1); ++i)
            {
               temp = abs(ABj[i]);
               if (temp > value)
               {
                  value = temp;
               }
               else if (isnan(temp))
                  return temp;
            }
            ABj += ldAB;
         }
      }
      else if ( normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t sum;
         complex<real_t> * ABj = AB;
         for (int_t j = 0; j < n; ++j)
         {
            sum = zero;
            for (int_t i = max(kU-j, intzero); i < min(n+kU-j, kL+kU+1); ++i)
            {
               sum += abs(ABj[i]);
            }
            if (sum > value)
            {
               value = sum;
            }
            else if (isnan(sum))
               return sum;
            ABj += ldAB;
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         real_t * Work = new real_t[n];
         real_t temp;
         for (int_t i = 0; i < n; ++i)
         {
            Work[i] = zero;
         }
         complex<real_t> * ABj = AB;
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = max(kU-j, intzero); i < min(n+kU-j, kL+kU+1); ++i)
            {
               Work[j+i-kU] += abs(ABj[i]);
            }
            ABj += ldAB;
         }
         for (int_t i = 0; i < n; ++i)
         {
            temp = Work[i];
            if (temp > value)
            {
               value = Work[i];
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
         complex<real_t> * ABj = AB;
         int_t L, K;
         for (int_t j = 0; j < n; ++j)
         {
            L = max(intzero,j-kU);
            K = kU -j+L;
            latl::lassq(min(n, j+kL+1)-L, ABj+K, 1, scale, sum);
            ABj += ldAB;
         }
         value = scale*sqrt(sum);
      }
      return value;
   }
}
#endif
