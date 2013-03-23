//
//  lantr.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/23/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lantr_h
#define _lantr_h

/// @file lantr.h  Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a triangular or trapezoidal matrix A.

#include "lassq.h"
#include "latl.h"
#include <cmath>

namespace LATL
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a real triangular or trapezoidal matrix A.
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
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper trapezoidal or lower trapezoidal.  The other part of A is not referenced.
   /// @param diag specifies whether or not A is unit triangular.  If diag = 'U' or 'u' then A is unit triangular; otherwise A is assumed to not be unit triangular.
   /// @param m Number of rows to be included in the norm. If m = n, then A is triangular.  m >= 0
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param A Real matrix size m-by-n.  
   /// @param ldA Column length of the matrix A.  ldA >= m
   /// @param Work Workspace vector of length m (optional).  If not used, workspace will be allocated and
   /// deallocated internally; only used for the infinity norm.
   /// @ingroup NORM
   
   template< typename real_t>
   real_t lantr(const char normType, const char uplo, const char diag, const int_t m, const int_t n, real_t * const A, const int_t ldA, real_t *Work=NULL)
   {
      using std::min;
      using std::max;
      using std::isnan;
      using std::abs;
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t value(0.0);
      if (m <= 0 || n <= 0)
         return value;

      if (normType == 'M' || normType == 'm')
      {
         real_t temp(0.0), * Aj = A;
         if (diag == 'U' || diag == 'u')
         {
            value = one;
            if (uplo == 'U' || uplo == 'u')
            {
               Aj += ldA;
               for (int_t j = 1; j < n; ++j)
               {
                  for (int_t i = 0; i < min(m, j); ++i)
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
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = j+1; i < m; ++i)
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
         else
         {
            value = zero;
            if (uplo == 'U' || uplo == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < min(m, j+1); ++i)
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
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = j; i < m; ++i)
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
      }
      else if (normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t sum(0.0), * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               if ((diag == 'U' || diag == 'u') && (j < m))
               {
                  sum = one;
                  for (int_t i = 0; i < j-1; ++i)
                  {
                     sum += abs(Aj[i]);
                  }
               }
               else
               {
                  sum = zero;
                  for (int_t i = 0; i < min(m, j+1); ++i)
                  {
                     sum += abs(Aj[i]);
                  }
               }
               if (value < sum)
                  value = sum;
               else if (isnan(sum))
                  return sum;
               Aj += ldA;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (diag == 'U' || diag == 'u')
               {
                  sum = one;
                  for (int_t i = j+1; i < m; ++i)
                  {
                     sum += abs(Aj[i]);
                  }
               }
               else
               {
                  sum = zero;
                  for (int_t i = j; i < m; ++i)
                  {
                     sum += abs(Aj[i]);
                  }
               }
               if (value < sum)
                  value = sum;
               else if (isnan(sum))
                  return sum;
               Aj += ldA;
            }
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         bool allocate=(Work==NULL)?1:0;
         if(allocate)
            Work = new real_t[m];
         real_t * Aj = A;
         real_t sum;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t i = 0; i < m; ++i)
               {
                  Work[i] = one;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < min(m, j); ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
                  Aj += ldA;
               }
            }
            else
            {
               for (int_t i = 0; i < m; ++i)
               {
                  Work[i] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < min(m, j+1); ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
                  Aj += ldA;
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = one;
               }
               for (int_t i = n; i < m; ++i)
                  Work[i] = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = j+1; i < m; ++i)
                     Work[i] += abs(Aj[i]);
                  Aj += ldA;
               }
            }
            else
            {
               for (int_t i = 0; i < m; ++i)
               {
                  Work[i] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = j; i < m; ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
                  Aj += ldA;
               }
            }
         }
         for (int_t i = 0; i < m; ++i)
         {
            sum = Work[i];
            if (value < sum)
               value = sum;
            else if (isnan(sum))
            {
               if(allocate)
                  delete [] Work;
               return sum;
            }
         }
         if(allocate)
            delete [] Work;
      }
      else if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale(1.0);
         real_t sum(1.0);
         real_t * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               sum = min(m, n);
               Aj += ldA;
               for (int_t j = 1; j < n; ++j)
               {
                  LATL::lassq(min(m, j), Aj, 1, scale, sum);
                  Aj += ldA;
               }
            }
            else
            {
               scale = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::lassq(min(m, j+1), Aj, 1, scale, sum);
                  Aj += ldA;
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               sum = min(m, n);
               Aj += 1;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::lassq(m-j-1, Aj, 1, scale, sum);
                  Aj += ldA+1;
               }
            }
            else
            {
               scale = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::lassq(m-j, Aj, 1, scale, sum);
                  Aj += ldA+1;
               }
            }
         }
         value = scale*sqrt(sum);
      }
      return value;
   }
   
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a complex triangular or trapezoidal matrix A.
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
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper trapezoidal or lower trapezoidal.  The other part of A is not referenced.
   /// @param diag specifies whether or not A is unit triangular.  If diag = 'U' or 'u' then A is unit triangular; otherwise A is assumed to not be unit triangular.
   /// @param m Number of rows to be included in the norm. If m = n, then A is triangular.  m >= 0
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param A Complex matrix size m-by-n.
   /// @param ldA Column length of the matrix A.  ldA >= m
   /// @param Work Workspace vector of length m (optional).  If not used, workspace will be allocated and
   /// deallocated internally; only used for the infinity norm.
   /// @ingroup NORM
   
   template< typename real_t>
   real_t lantr(const char normType, const char uplo, const char diag, const int_t m, const int_t n, complex<real_t> * const A, const int_t ldA, real_t *Work=NULL)
   {
      using std::min;
      using std::max;
      using std::isnan;
      using std::abs;
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t value(0.0);
      if (m <= 0 || n <= 0)
         return value;

      if (normType == 'M' || normType == 'm')
      {
         real_t temp(0.0);
         complex<real_t> * Aj = A;
         if (diag == 'U' || diag == 'u')
         {
            value = one;
            if (uplo == 'U' || uplo == 'u')
            {
               Aj += ldA;
               for (int_t j = 1; j < n; ++j)
               {
                  for (int_t i = 0; i < min(m, j); ++i)
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
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = j+1; i < m; ++i)
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
         else
         {
            value = zero;
            if (uplo == 'U' || uplo == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < min(m, j+1); ++i)
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
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = j; i < m; ++i)
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
      }
      else if (normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t sum(0.0);
         complex<real_t> * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               if ((diag == 'U' || diag == 'u') && (j < m))
               {
                  sum = one;
                  for (int_t i = 0; i < j-1; ++i)
                  {
                     sum += abs(Aj[i]);
                  }
               }
               else
               {
                  sum = zero;
                  for (int_t i = 0; i < min(m, j+1); ++i)
                  {
                     sum += abs(Aj[i]);
                  }
               }
               if (value < sum)
                  value = sum;
               else if (isnan(sum))
                  return sum;
               Aj += ldA;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (diag == 'U' || diag == 'u')
               {
                  sum = one;
                  for (int_t i = j+1; i < m; ++i)
                  {
                     sum += abs(Aj[i]);
                  }
               }
               else
               {
                  sum = zero;
                  for (int_t i = j; i < m; ++i)
                  {
                     sum += abs(Aj[i]);
                  }
               }
               if (value < sum)
                  value = sum;
               else if (isnan(sum))
                  return sum;
               Aj += ldA;
            }
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         bool allocate=(Work==NULL)?1:0;
         if(allocate)
            Work = new real_t[m];
         complex<real_t> * Aj = A;
         real_t sum;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t i = 0; i < m; ++i)
               {
                  Work[i] = one;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < min(m, j); ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
                  Aj += ldA;
               }
            }
            else
            {
               for (int_t i = 0; i < m; ++i)
               {
                  Work[i] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < min(m, j+1); ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
                  Aj += ldA;
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = one;
               }
               for (int_t i = n; i < m; ++i)
                  Work[i] = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = j+1; i < m; ++i)
                     Work[i] += abs(Aj[i]);
                  Aj += ldA;
               }
            }
            else
            {
               for (int_t i = 0; i < m; ++i)
               {
                  Work[i] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = j; i < m; ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
                  Aj += ldA;
               }
            }
         }
         for (int_t i = 0; i < m; ++i)
         {
            sum = Work[i];
            if (value < sum)
               value = sum;
            else if (isnan(sum))
            {
               if(allocate)
                  delete [] Work;
               return sum;
            }
         }
         if(allocate)
            delete [] Work;
      }
      else if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale(1.0);
         real_t sum(1.0);
         complex<real_t> * Aj = A;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               sum = min(m, n);
               Aj += ldA;
               for (int_t j = 1; j < n; ++j)
               {
                  LATL::lassq(min(m, j), Aj, 1, scale, sum);
                  Aj += ldA;
               }
            }
            else
            {
               scale = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::lassq(min(m, j+1), Aj, 1, scale, sum);
                  Aj += ldA;
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               sum = min(m, n);
               Aj += 1;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::lassq(m-j-1, Aj, 1, scale, sum);
                  Aj += ldA+1;
               }
            }
            else
            {
               scale = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::lassq(m-j, Aj, 1, scale, sum);
                  Aj += ldA+1;
               }
            }
         }
         value = scale*sqrt(sum);
      }
      return value;
   }
}
#endif
