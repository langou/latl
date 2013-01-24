//
//  lanhb.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/23/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lanhb_h
#define _lanhb_h

/// @file lanhb.h Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a Hermitian band matrix A.

#include "lassq.h"
#include "latl.h"
#include <cmath>

namespace latl
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a Hermitian band matrix.
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
   /// @param uplo Indicates whether the symmetric band matrix A is stored as upper triangular or lower triangular.
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param k The number of super- or sub-diagonals within the band of A.  k >= 0.
   /// @param AB Complex matrix size ldAB-by-n.  On entry, the matrix A in band storage.
   /// @param ldAB Column length of the matrix AB.  ldAB >= k+1
   /// @ingroup NORM
   
   
   template< typename real_t>
   real_t lanhb(const char normType, const char uplo, const int_t n, const int_t k, complex<real_t> * const AB, const int_t ldAB)
   {
      using std::isnan;
      using std::abs;
      using std::min;
      using std::max;
      using std::real;
      const real_t zero(0.0);
      const int_t intzero(0);
      real_t value(0.0);
      if (n == 0)
         return value;
      if (normType == 'M' || normType == 'm')
      {
         complex<real_t> * ABj = AB;
         real_t temp(0.0);
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               for (int_t i = max(k-j, intzero); i < k; ++i)
               {
                  temp = abs(ABj[i]);
                  if ((temp > value) || (isnan(temp)))
                     value = temp;
               }
               temp = abs(real(ABj[k]));
               if (value < temp || isnan(temp))
                  value = temp;
               ABj += ldAB;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               temp = abs(real(ABj[0]));
               if ((temp > value) || (isnan(temp)))
                  value = temp;
               for (int_t i = 1; i < min(n-j, k+1); ++i)
               {
                  temp = abs(ABj[i]);
                  if ((temp > value) || (isnan(temp)))
                  {
                     value = temp;
                  }
               }
               ABj += ldAB;
            }
         }
      }
      else if ( normType == 'O' || normType == 'o' || normType == '1' || normType == 'I' || normType == 'i')
      {
         real_t * Work = new real_t[n];
         complex<real_t> * ABj = AB;
         real_t sum(0.0);
         for (int_t i = 0; i < n; ++i)
         {
            Work[i] = zero;
         }
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               sum = zero;
               for (int_t i = max(k-j, intzero); i < k; ++i)
               {
                  sum += abs(ABj[i]);
                  Work[j+i-k] += abs(ABj[i]);
               }
               Work[j] = sum + abs(real(ABj[k]));
               ABj += ldAB;
            }
            for (int_t i = 0; i < n; ++i)
            {
               sum = Work[i];
               if (sum > value || isnan(sum))
                  value = sum;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               sum = Work[j] + abs(real(ABj[0]));
               for (int_t i = 1; i < min(n-j, k+1); ++i)
               {
                  sum += abs(ABj[i]);
                  Work[j+i] += abs(ABj[i]);
               }
               if(value < sum || isnan(sum))
                  value = sum;
               ABj += ldAB;
            }
         }
         delete [] Work;
      }
      else if ( normType == 'F' || normType == 'f' || normType == 'E' ||normType == 'e')
      {
         real_t scale(0.0), sum(1.0), temp(0.0);
         complex<real_t> * ABj = AB;
         int_t L = 0;
         if (k > 0)
         {
            if (uplo == 'U' || uplo == 'u')
            {
               ABj += ldAB + k-1;
               for (int_t j = 1; j < n; ++j)
               {
                  latl::lassq(min(j, k), ABj, 1, scale, sum);
                  if (k-j > 0)
                     ABj += ldAB-1;
                  else
                     ABj += ldAB;
               }
               L = k;
            }
            else
            {
               ABj += 1;
               for (int_t j = 0; j < n-1; ++j)
               {
                  latl::lassq(min(k, n-j-1), ABj, 1, scale, sum);
                  ABj += ldAB;
               }
            }
            sum *= 2;
         }
         ABj = AB;
         for (int_t j = 0; j < n; ++j)
         {
            if (real(ABj[L]) != zero)
            {
               temp = abs(real(ABj[L]));
               if (scale < temp)
               {
                  sum = 1.0 + sum*((scale/temp)*(scale/temp));
                  scale = temp;
               }
               else
               {
                  sum += (temp/scale)*(temp/scale);
               }
            }
            ABj += ldAB;
         }
         value = scale*sqrt(sum);
      }
      return value;
   }
}
#endif
