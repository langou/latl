//
//  lanhe.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/24/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lanhe_h
#define _lanhe_h

/// @file lanhe.h Returns the norm of a Hermitian matrix.

#include "lassq.h"
#include "latl.h"
#include <cmath>

namespace LATL
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest
   /// absolute value of a complex Hermitian matrix A.
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
   /// @param uplo Indicates whether the Hermitian matrix A is stored as upper triangular or lower triangular.
   /// The other part of A is not referenced.
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param A Complex matrix size n-by-n.
   /// @param ldA Column length of the matrix A.  ldA >= n
   /// @ingroup AUX
   
   template<typename real_t>
   real_t LANHE(const char normType, const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA)
   {
      using std::abs;
      using std::real;
      using std::isnan;
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t value(0.0);
      if (n == 0)
         return value;

      if (normType == 'M' || normType == 'm')
      {
         real_t temp(0.0);
         complex<real_t> * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               for (int_t i = 0; i < j; ++i)
               {
                  temp = abs(Aj[i]);
                  if (value < temp)
                     value = temp;
                  else if (isnan(temp))
                     return temp;
               }
               temp = abs(real(Aj[j]));
               if (value < temp)
                  value = temp;
               else if (isnan(temp))
                  return temp;
               Aj += ldA;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               temp = abs(real(Aj[j]));
               if (value < temp)
                  value = temp;
               else if (isnan(temp))
                  return temp;
               for (int_t i = j+1; i < n; ++i)
               {
                  temp = abs(Aj[i]);
                  if (value < temp)
                     value = temp;
                  else if (isnan(temp))
                     return temp;
               }
               Aj += ldA;
            }
         }
      }
      else if (normType == 'O' ||normType == 'o' || normType == 'I' || normType == 'i' || normType == '1')
      {
         real_t sum(0.0), temp(0.0);
         complex<real_t> * Aj = A;
         real_t * Work = new real_t[n];
         for (int_t i = 0; i < n; ++i)
            Work[i] = zero;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               sum = zero;
               for (int_t i = 0; i < j; ++i)
               {
                  temp = abs(Aj[i]);
                  sum += temp;
                  Work[i] += temp;
               }
               Work[j] = sum + abs(real(Aj[j]));
               Aj += ldA;
            }
            for (int_t i = 0; i < n; ++i)
            {
               sum = Work[i];
               if (value < sum)
                  value = sum;
               else if (isnan(sum))
               {
                  delete [] Work;
                  return sum;
               }
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               sum = Work[j] + abs(real(Aj[j]));
               for (int_t i = j+1; i < n; ++i)
               {
                  temp = abs(Aj[i]);
                  sum += temp;
                  Work[i] += temp;
               }
               Aj += ldA;
            }
            if (value < sum)
               value = sum;
            else if (isnan(sum))
            {
               delete [] Work;
               return sum;
            }
         }
         delete [] Work;
      }
      else if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale(0.0);
         real_t sum(1.0);
         real_t temp(0.0), temp2(0.0);
         complex<real_t> * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            Aj += ldA;
            for (int_t j = 1; j < n; ++j)
            {
               LATL::LASSQ(j, Aj, 1, scale, sum);
               Aj += ldA;
            }
         }
         else
         {
            Aj += 1;
            for (int_t j = 0; j < n-1; ++j)
            {
               LATL::LASSQ(n-j-1, Aj, 1, scale, sum);
               Aj += ldA+1;
            }
         }
         sum *= 2;
         Aj = A;
         for (int_t i = 0; i < n; ++i)
         {
            if (real(Aj[i]) != zero)
            {
               temp = abs(real(Aj[i]));
               if (scale < temp)
               {
                  temp2 = (scale/temp);
                  sum = one + sum*(temp2*temp2);
                  scale = temp;
               }
               else
               {
                  temp2 = temp/scale;
                  sum += (temp2*temp2);
               }
            }
            Aj += ldA;
         }
         value = scale*sqrt(sum);
      }
      return value;
   }
}

#endif
