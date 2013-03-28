//
//  lansy.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/20/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lansy_h
#define _lansy_h

/// @file lansy.h Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest
/// absolute value of a symmetric matrix A.

#include "latl.h"
#include <cmath>
#include "lassq.h"

namespace LATL
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest
   /// absolute value of a real symmetric matrix.
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
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.
   /// The other triangular part of A is not referenced.
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param A Real symmetric matrix size ldA-by-n.
   /// @param ldA Column length of the matrix A.  ldA >= n
   /// @ingroup AUX
   
   template< typename real_t>
   real_t LANSY(const char normType, const char uplo, const int_t n, real_t * const A, const int_t ldA)
   {
      using std::isnan;
      using std::abs;
      using std::sqrt;
      const real_t zero(0.0);
      real_t value(0.0);
      if (n == 0)
         return zero;

      if (normType == 'M' || normType == 'm')
      {
         real_t * Aj = A;
         real_t temp = 0;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               for (int_t i = 0; i <= j; ++i)
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
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               for (int_t i = j; i < n; ++i)
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
      }
      else if ( normType == 'O' || normType == 'o' || normType == '1' || normType == 'I' || normType == 'i')
      {
         real_t sum, * Aj = A;
         real_t * Work = new real_t[n];
         for (int_t i = 0; i < n; ++i)
         {
            Work[i] = zero;
         }
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               sum = zero;
               for (int_t i = 0; i < j; ++i)
               {
                  sum += abs(Aj[i]);
                  Work[i] += abs(Aj[i]);
               }
               Work[j] = sum + abs(Aj[j]);
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
               sum = Work[j] + abs(Aj[j]);
               for (int_t i = j+1; i < n; ++i)
               {
                  sum += abs(Aj[i]);
                  Work[i] += abs(Aj[i]);
               }
               Aj += ldA;
               if (value < sum)
                  value = sum;
               else if (isnan(sum))
               {
                  delete [] Work;
                  return sum;
               }
            }
         }
         delete [] Work;
      }
      else if ( normType == 'F' || normType == 'f' || normType == 'E' ||normType == 'e')
      {
         real_t scale(0.0), sum(1.0), * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 1; j < n; ++j)
            {
               Aj += ldA;
               LATL::LASSQ(j, Aj, 1, scale, sum);
            }
         }
         else
         {
            for (int_t j = 0; j < n-1; ++j)
            {
               Aj += 1;
               LATL::LASSQ(n-j-1, Aj, 1, scale, sum);
               Aj += ldA;
            }
         }
         sum *= 2;
         LATL::LASSQ(n, A, ldA+1, scale, sum);
         value = scale*sqrt(sum);
      }
      return value;
   }

   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute
   /// value of a complex symmetric matrix.
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
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.
   /// The other triangular part of A is not referenced.
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param A Complex symmetric matrix size ldA-by-n.
   /// @param ldA Column length of the matrix A.  ldA >= n
   /// @ingroup AUX
   
   template< typename real_t>
   real_t LANSY(const char normType, const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA)
   {
      using std::isnan;
      using std::abs;
      using std::sqrt;
      const real_t zero(0.0);
      real_t value(0.0);
      if (n == 0)
         return zero;

      if (normType == 'M' || normType == 'm')
      {
         complex<real_t> * Aj = A;
         real_t temp = 0;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               for (int_t i = 0; i <= j; ++i)
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
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               for (int_t i = j; i < n; ++i)
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
      }
      else if ( normType == 'O' || normType == 'o' || normType == '1' || normType == 'I' || normType == 'i')
      {
         real_t sum;
         complex<real_t> * Aj = A;
         real_t * Work = new real_t[n];
         for (int_t i = 0; i < n; ++i)
         {
            Work[i] = zero;
         }
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               sum = zero;
               for (int_t i = 0; i < j; ++i)
               {
                  sum += abs(Aj[i]);
                  Work[i] += abs(Aj[i]);
               }
               Work[j] = sum + abs(Aj[j]);
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
               sum = Work[j] + abs(Aj[j]);
               for (int_t i = j+1; i < n; ++i)
               {
                  sum += abs(Aj[i]);
                  Work[i] += abs(Aj[i]);
               }
               Aj += ldA;
               if (value < sum)
                  value = sum;
               else if (isnan(sum))
               {
                  delete [] Work;
                  return sum;
               }
            }
         }
         delete [] Work;
      }
      else if ( normType == 'F' || normType == 'f' || normType == 'E' ||normType == 'e')
      {
         real_t scale(0.0), sum(1.0);
         complex<real_t> * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 1; j < n; ++j)
            {
               Aj += ldA;
               LATL::LASSQ(j, Aj, 1, scale, sum);
            }
         }
         else
         {
            for (int_t j = 0; j < n-1; ++j)
            {
               Aj += 1;
               LATL::LASSQ(n-j-1, Aj, 1, scale, sum);
               Aj += ldA;
            }
         }
         sum *= 2;
         LATL::LASSQ(n, A, ldA+1, scale, sum);
         value = scale*sqrt(sum);
      }
      return value;
   }
}


#endif
