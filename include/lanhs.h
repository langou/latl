//
//  lanhs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/20/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lanhs_h
#define _lanhs_h

/// @file lanhs.h Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of an upper Hessenberg matrix A.

#include "latl.h"
#include "lassq.h"

namespace LATL
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a real upper Hessenberg matrix A.
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
   /// @param A Real matrix size ldA-by-n.  Elements of A below the first subdiagonal are not referenced.
   /// @param ldA Column length of the matrix A.  ldA >= n
   /// @param Work Workspace vector of length m (optional).  If not used, workspace will be allocated and
   /// deallocated internally; only used for the infinity norm.
   /// @ingroup NORM
   
   template <typename real_t>
   real_t lanhs( const char normType, const int_t n, real_t * A, const int_t ldA, real_t *Work=NULL)
   {
      using std::min;
      using std::isnan;
      using std::abs;
      real_t value(0.0);
      const real_t zero(0.0);
      if (n <= 0)
         return value;

      if (normType == 'M' || normType == 'm')
      {
         real_t temp(0.0), * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = 0; i < min(n, j+2); ++i)
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
      else if (normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t sum(0.0), * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            sum = zero;
            for (int_t i = 0; i < min(n, j+2); ++i)
            {
               sum += abs(Aj[i]);
            }
            if (value < sum)
               value = sum;
            else if (isnan(sum))
               return sum;
            Aj += ldA;
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         if (n == 1)
            value = abs(A[0]);
         else
         {
            bool allocate=(Work==NULL)?1:0;
            if(allocate)
               Work = new real_t[n];
            real_t * Aj = A;
            real_t temp;
            for (int_t i = 0; i < n; ++i)
            {
               Work[i] = zero;
            }
            for (int_t j = 0; j < n; ++j)
            {
               for (int_t i = 0; i < min(n, j+2); ++i)
               {
                  Work[i] += abs(Aj[i]);
               }
               Aj += ldA;
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
      }
      else if (normType == 'E' || normType == 'F' || normType == 'e' || normType == 'f')
      {
         real_t scale(0.0);
         real_t sum(1.0);
         real_t * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            LATL::lassq(min(n, j+2), Aj, 1, scale, sum);
            Aj += ldA;
         }
         value = scale * sqrt(sum);
      }
      return value;
   }
   
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a complex upper Hessenberg matrix A.
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
   /// @param A Complex matrix size ldA-by-n.  Elements of A below the first subdiagonal are not referenced.
   /// @param ldA Column length of the matrix A.  ldA >= n
   /// @param Work Workspace vector of length m (optional).  If not used, workspace will be allocated and
   /// deallocated internally; only used for the infinity norm.
   /// @ingroup NORM
   
   template <typename real_t>
   real_t lanhs( const char normType, const int_t n, complex<real_t> * A, const int_t ldA, real_t *Work=NULL)
   {
      using std::min;
      using std::isnan;
      using std::abs;
      real_t value(0.0);
      const real_t zero(0.0);
      if (n <= 0)
         return value;

      if (normType == 'M' || normType == 'm')
      {
         real_t temp(0.0);
         complex<real_t> * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = 0; i < min(n, j+2); ++i)
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
      else if (normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t sum(0.0);
         complex<real_t> * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            sum = zero;
            for (int_t i = 0; i < min(n, j+2); ++i)
            {
               sum += abs(Aj[i]);
            }
            if (value < sum)
               value = sum;
            else if (isnan(sum))
               return sum;
            Aj += ldA;
         }
      }
      else if (normType == 'I' || normType == 'i')
      {
         if (n == 1)
            value = abs(A[0]);
         else
         {
            bool allocate=(Work==NULL)?1:0;
            if(allocate)
               Work = new real_t[n];
            complex<real_t> * Aj = A;
            real_t temp;
            for (int_t i = 0; i < n; ++i)
            {
               Work[i] = zero;
            }
            for (int_t j = 0; j < n; ++j)
            {
               for (int_t i = 0; i < min(n, j+2); ++i)
               {
                  Work[i] += abs(Aj[i]);
               }
               Aj += ldA;
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
      }
      else if (normType == 'E' || normType == 'F' || normType == 'e' || normType == 'f')
      {
         real_t scale(0.0);
         real_t sum(1.0);
         complex<real_t> * Aj = A;
         for (int_t j = 0; j < n; ++j)
         {
            LATL::lassq(min(n, j+2), Aj, 1, scale, sum);
            Aj += ldA;
         }
         value = scale * sqrt(sum);
      }
      return value;
   }
   
}

#endif
