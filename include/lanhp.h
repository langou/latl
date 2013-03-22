//
//  lanhp.h
//  Linear Algebra Template Library  
//
//  Created by Stephanie Patterson on 1/26/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lanhp_h
#define _lanhp_h

/// @file lanhp.h Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a general matrix A.

#include "lassq.h"
#include <cmath>
#include "latl.h"

namespace latl
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a complex Hermitian matrix A supplied in packed form.
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
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param AP Pointer to packed complex Hermitian n-by-n matrix A.
   /// @param Work Workspace vector of length n (optional).  If not used, workspace will be allocated and
   /// deallocated internally; only used for the infinity norm.

   /// @ingroup NORM
   
   template< typename real_t>
   real_t lanhp(const char normType, const char uplo, const int_t n, complex<real_t> * const AP, real_t *Work=NULL)
   {
      using std::abs;
      using std::isnan;
      const real_t zero(0.0);
      const real_t one(1.0);
      real_t value(0.0);
      if (n == 0)
         return value;
      bool allocate=(Work==NULL)?1:0;

      if (normType == 'M' || normType == 'm')
      {
         real_t temp;
         complex<real_t> * Aj = AP;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
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
               temp = abs(real(Aj[j]));
               if (value < temp)
                  value = temp;
               else if (isnan(temp))
                  return temp;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               temp = abs(real(Aj[0]));
               if (value < temp)
                  value = temp;
               else if (isnan(temp))
                  return temp;
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
      else if (normType == 'O' || normType == 'o' || normType == '1' || normType == 'I' || normType == 'i')
      {
         complex<real_t> * Aj = AP;
         if(allocate)
            Work = new real_t[n];
         real_t temp, sum;
         for (int_t i = 0; i < n; ++i)
         {
            Work[i] = zero;
         }
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               sum = zero;
               Aj += j;
               for (int_t i = 0; i < j; ++i)
               {
                  temp = abs(Aj[i]);
                  Work[i] += temp;
                  sum += temp;
               }
               Work[j] = sum + abs(real(Aj[j]));
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
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               sum = Work[j] + abs(real(Aj[0]));
               for (int_t i = 1; i < n-j; ++i)
               {
                  temp = abs(Aj[i]);
                  sum += temp;
                  Work[i+j] += temp;
               }
               if (value < sum)
                  value = sum;
               else if (isnan(sum))
               {
                  if(allocate)
                     delete [] Work;
                  return sum;
               }
               Aj += n-j;
            }
         }
         if(allocate)
            delete [] Work;
      }
      else if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale(0.0);
         real_t sum(1.0);
         complex<real_t> * Aj = AP;
         real_t temp, temp2;
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 1; j < n; ++j)
            {
               Aj += j;
               latl::lassq(j, Aj, 1, scale, sum);
            }
            sum *= 2;
            Aj = AP;
            for (int_t j = 0; j < n; ++j)
            {
               Aj += j;
               if (real(Aj[j]) != zero)
               {
                  temp = abs(real(Aj[j]));
                  if (scale < temp)
                  {
                     temp2 = scale/temp;
                     sum = one + sum*(temp2*temp2);
                     scale = temp;
                  }
                  else
                  {
                     temp2 = temp/scale;
                     sum += (temp2*temp2);
                  }
               }
            }
         }
         else
         {
            Aj += 1;
            for (int_t j = 0; j < n-1; ++j)
            {
               latl::lassq(n-j-1, Aj, 1, scale, sum);
               Aj += n-j;
            }
            sum *= 2;
            Aj = AP;
            for (int_t j = 0; j < n; ++j)
            {
               if (real(Aj[0]) != zero)
               {
                  temp = abs(real(Aj[0]));
                  if (scale < temp)
                  {
                     temp2 = scale/temp;
                     sum = one + sum*(temp2*temp2);
                     scale = temp;
                  }
                  else
                  {
                     temp2 = temp/scale;
                     sum += (temp2*temp2);
                  }
               }
               Aj += n-j;
            }
         }
         value = scale*sqrt(sum);
      }
      return value;
   }
}

#endif
