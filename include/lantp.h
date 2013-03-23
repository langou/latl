//
//  lantp.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/25/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lantp_h
#define _lantp_h

/// @file lantp.h  Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a triangular matrix A.

#include "lassq.h"
#include <cmath>
#include "latl.h"

namespace latl
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a real triangular matrix A supplied in packed form.
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
   /// @param uplo Specifies whether A is stored as upper or lower triangular.  If uplo = 'U' or 'u' then A is upper triangular; otherwise A is assumed to be lower triangular.
   /// @param diag specifies whether or not A is unit triangular.  If diag = 'U' or 'u' then A is unit triangular and the diagonal elements of A are not referenced; otherwise A is assumed to not be unit triangular.
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param AP Pointer to packed real triangular n-by-n matrix A.
   /// @param Work Workspace vector of length m (optional).  If not used, workspace will be allocated and
   /// deallocated internally; only used for the infinity norm.
   /// @ingroup NORM

   template< typename real_t>
   real_t lantp(const char normType, const char uplo, const char diag, const int_t n, real_t * const AP, real_t *Work=NULL)
   {
      using std::abs;
      using std::isnan;
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t value(0.0);
      if (n == 0)
         return value;

      if (normType == 'M' || normType == 'm')
      {
         real_t temp;
         real_t * Aj = AP;
         if (diag == 'U' || diag == 'u')
         {
            value = one;
            if (uplo == 'U' || uplo == 'u')
            {
               for (int_t j = 1; j < n; ++j)
               {
                  Aj += j;
                  for (int_t i = 0; i < j; ++i)
                  {
                     temp = abs(Aj[i]);
                     if (value < temp)
                        value = temp;
                     else if (isnan(temp))
                        return temp;
                  }
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 1; i < n-j; ++i)
                  {
                     temp = abs(Aj[i]);
                     if (value < temp)
                        value = temp;
                     else if (isnan(temp))
                        return temp;
                  }
                  Aj += n-j;
               }
            }
         }
         else
         {
            if (uplo == 'U' || uplo == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  for (int_t i = 0; i <= j; ++i)
                  {
                     temp = abs(Aj[i]);
                     if (value < temp)
                        value = temp;
                     else if (isnan(temp))
                        return temp;
                  }
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < n-j; ++i)
                  {
                     temp = abs(Aj[i]);
                     if (value < temp)
                        value = temp;
                     else if (isnan(temp))
                        return temp;
                  }
                  Aj += n-j;
               }
            }
         }
      }
      else if (normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t * Aj = AP;
         real_t sum;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  sum = one;
                  for (int_t i = 0; i < j; ++i)
                     sum += abs(Aj[i]);
                  if (value < sum)
                     value = sum;
                  else if (isnan(sum))
                     return sum;
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  sum = zero;
                  for (int_t i = 0; i <= j; ++i)
                     sum += abs(Aj[i]);
                  if (value < sum)
                     value = sum;
                  else if (isnan(sum))
                     return sum;
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  sum = one;
                  for (int_t i = 1; i < n-j; ++i)
                     sum += abs(Aj[i]);
                  if (value < sum)
                     value = sum;
                  else if (isnan(sum))
                     return sum;
                  Aj += n-j;
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  sum = zero;
                  for (int_t i = 0; i < n-j; ++i)
                     sum += abs(Aj[i]);
                  if (value < sum)
                     value = sum;
                  else if (isnan(sum))
                     return sum;
                  Aj += n-j;
               }
            }
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         real_t * Aj = AP;
         bool allocate=(Work==NULL)?1:0;
         if(allocate)
            Work = new real_t[n];
         real_t temp;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = one;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  for (int_t i = 0; i < j; ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
               }
            }
            else
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  for (int_t i = 0; i <= j; ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
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
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 1; i < n-j; ++i)
                  {
                     Work[i+j] += abs(Aj[i]);
                  }
                  Aj += n-j;
               }
            }
            else
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < n-j; ++i)
                  {
                     Work[i+j] += abs(Aj[i]);
                  }
                  Aj += n-j;
               }
            }
         }
         for (int_t i = 0; i < n; ++i)
         {
            temp = Work[i];
            if (value < temp)
               value = temp;
            else if (isnan(temp))
            {
               if(allocate)
                  delete [] Work;
               return temp;
            }
         }
         if(allocate)
            delete [] Work;
      }
      else if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale(1.0);
         real_t sum(1.0);
         real_t * Aj = AP;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               sum = n;
               for (int_t j = 1; j < n; ++j)
               {
                  Aj += j;
                  latl::lassq(j, Aj, 1, scale, sum);
               }
            }
            else
            {
               scale = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  latl::lassq(j+1, Aj, 1, scale, sum);
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               sum = n;
               Aj += 1;
               for (int_t j = 0; j < n-1; ++j)
               {
                  latl::lassq(n-j-1, Aj, 1, scale, sum);
                  Aj += n-j;
               }
            }
            else
            {
               scale = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  latl::lassq(n-j, Aj, 1, scale, sum);
                  Aj += n-j;
               }
            }
         }
         value = scale*sqrt(sum);
      }
      return value;      
   }
   
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a complex triangular matrix A supplied in packed form.
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
   /// @param uplo Specifies whether A is stored as upper or lower triangular.  If uplo = 'U' or 'u' then A is upper triangular; otherwise A is assumed to be lower triangular.
   /// @param diag specifies whether or not A is unit triangular.  If diag = 'U' or 'u' then A is unit triangular and the diagonal elements of A are not referenced; otherwise A is assumed to not be unit triangular.
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param AP Pointer to packed complex triangular n-by-n matrix A.
   /// @param Work Workspace vector of length m (optional).  If not used, workspace will be allocated and
   /// deallocated internally; only used for the infinity norm.
   /// @ingroup NORM
   
   template< typename real_t>
   real_t lantp(const char normType, const char uplo, const char diag, const int_t n, complex<real_t> * const AP, real_t *Work=NULL)
   {
      using std::abs;
      using std::isnan;
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t value(0.0);
      if (n == 0)
         return value;

      if (normType == 'M' || normType == 'm')
      {
         real_t temp;
         complex<real_t> * Aj = AP;
         if (diag == 'U' || diag == 'u')
         {
            value = one;
            if (uplo == 'U' || uplo == 'u')
            {
               for (int_t j = 1; j < n; ++j)
               {
                  Aj += j;
                  for (int_t i = 0; i < j; ++i)
                  {
                     temp = abs(Aj[i]);
                     if (value < temp)
                        value = temp;
                     else if (isnan(temp))
                        return temp;
                  }
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 1; i < n-j; ++i)
                  {
                     temp = abs(Aj[i]);
                     if (value < temp)
                        value = temp;
                     else if (isnan(temp))
                        return temp;
                  }
                  Aj += n-j;
               }
            }
         }
         else
         {
            if (uplo == 'U' || uplo == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  for (int_t i = 0; i <= j; ++i)
                  {
                     temp = abs(Aj[i]);
                     if (value < temp)
                        value = temp;
                     else if (isnan(temp))
                        return temp;
                  }
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < n-j; ++i)
                  {
                     temp = abs(Aj[i]);
                     if (value < temp)
                        value = temp;
                     else if (isnan(temp))
                        return temp;
                  }
                  Aj += n-j;
               }
            }
         }
      }
      else if (normType == 'O' || normType == 'o' || normType == '1')
      {
         complex<real_t> * Aj = AP;
         real_t sum;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  sum = one;
                  for (int_t i = 0; i < j; ++i)
                     sum += abs(Aj[i]);
                  if (value < sum)
                     value = sum;
                  else if (isnan(sum))
                     return sum;
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  sum = zero;
                  for (int_t i = 0; i <= j; ++i)
                     sum += abs(Aj[i]);
                  if (value < sum)
                     value = sum;
                  else if (isnan(sum))
                     return sum;
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  sum = one;
                  for (int_t i = 1; i < n-j; ++i)
                     sum += abs(Aj[i]);
                  if (value < sum)
                     value = sum;
                  else if (isnan(sum))
                     return sum;
                  Aj += n-j;
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  sum = zero;
                  for (int_t i = 0; i < n-j; ++i)
                     sum += abs(Aj[i]);
                  if (value < sum)
                     value = sum;
                  else if (isnan(sum))
                     return sum;
                  Aj += n-j;
               }
            }
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         complex<real_t> * Aj = AP;
         bool allocate=(Work==NULL)?1:0;
         if(allocate)
            Work = new real_t[n];
         real_t temp;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = one;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  for (int_t i = 0; i < j; ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
               }
            }
            else
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  for (int_t i = 0; i <= j; ++i)
                  {
                     Work[i] += abs(Aj[i]);
                  }
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
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 1; i < n-j; ++i)
                  {
                     Work[i+j] += abs(Aj[i]);
                  }
                  Aj += n-j;
               }
            }
            else
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < n-j; ++i)
                  {
                     Work[i+j] += abs(Aj[i]);
                  }
                  Aj += n-j;
               }
            }
         }
         for (int_t i = 0; i < n; ++i)
         {
            temp = Work[i];
            if (value < temp)
               value = temp;
            else if (isnan(temp))
            {
               if(allocate)
                  delete [] Work;
               return temp;
            }
         }
         if(allocate)
            delete [] Work;
      }
      else if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale(1.0);
         real_t sum(1.0);
         complex<real_t> * Aj = AP;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               sum = n;
               for (int_t j = 1; j < n; ++j)
               {
                  Aj += j;
                  latl::lassq(j, Aj, 1, scale, sum);
               }
            }
            else
            {
               scale = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  Aj += j;
                  latl::lassq(j+1, Aj, 1, scale, sum);
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               sum = n;
               Aj += 1;
               for (int_t j = 0; j < n-1; ++j)
               {
                  latl::lassq(n-j-1, Aj, 1, scale, sum);
                  Aj += n-j;
               }
            }
            else
            {
               scale = zero;
               for (int_t j = 0; j < n; ++j)
               {
                  latl::lassq(n-j, Aj, 1, scale, sum);
                  Aj += n-j;
               }
            }
         }
         value = scale*sqrt(sum);
      }
      return value;
   }
}

#endif
