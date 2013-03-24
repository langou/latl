//
//  langt.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/15/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _langt_h
#define _langt_h

/// @file langt.h Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a tridiagonal matrix A.

#include "lassq.h"
#include <cmath>
#include "latl.h"

namespace LATL
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a real tridiagonal matrix.
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
   /// @param n Number of columns to be included in the norm.  n >= 0
   /// @param DL Real array, size n-1.  On entry, DL must contain the (n-1) subdiagonal elements of A.
   /// @param D Real array, size n.  On entry, D must contain the diagonal elements of A.
   /// @param DU Real array, size (n-1).  On entry, DU must contain the (n-1) superdiagonal elements of A.
   /// @ingroup NORM

   template <typename real_t>
   real_t LANGT(const char normType, const int_t n, real_t * const DL, real_t * const D, real_t * const DU)
   {
      using std::isnan;
      using std::abs;
      using std::sqrt;
      const real_t zero = 0.0;
      const real_t one = 1.0;
      real_t value = zero;
      if (n <= 0)
         return value;
      if (normType == 'M' || normType == 'm')
      {
         real_t temp;
         value = abs(D[n-1]);
         for (int_t i = 0; i < n-1; ++i)
         {
            temp = abs(DL[i]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
            temp = abs(D[i]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
            temp = abs(DU[i]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
         }
      }
      if (normType == 'O' || normType == 'o' || normType == '1')
      {
         if (n == 1)
            value = abs(D[0]);
         else
         {
            real_t temp;
            
            value = abs(D[0]) + abs(DL[0]);
            temp = abs(D[n-1]) + abs(DU[n-2]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
            for (int_t i = 1; i < n-1; ++i)
            {
               temp = abs(D[i])+abs(DU[i-1])+abs(DL[i]);
               if (value < temp)
                  value = temp;
               else if (isnan(temp))
                  return temp;
            }
         }
      }
      if (normType == 'I' || normType == 'i')
      {
         if (n == 1)
            value = abs(D[0]);
         else
         {
            real_t temp;
            
            value = abs(D[0]) + abs(DU[0]);
            temp = abs(D[n-1]) + abs(DL[n-2]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
            for (int_t i = 1; i < n-1; ++i)
            {
               temp = abs(D[i]) + abs(DL[i-1]) + abs(DU[i]);
               if (value < temp)
                  value = temp;
               else if (isnan(temp))
                  return temp;
            }
         }
      }
      if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale = zero;
         real_t sum = one;
         
         LATL::LASSQ(n, D, 1, scale, sum);
         if (n > 1)
         {
            LATL::LASSQ(n-1, DL, 1, scale, sum);
            LATL::LASSQ(n-1, DU, 1, scale, sum);
         }
         value = scale*sqrt(sum);
      }
      return value;
   }

   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a complex tridiagonal matrix.
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
   /// @param n Number of columns to be included in the norm.  n >= 0
   /// @param DL Real array, size n-1.  On entry, DL must contain the (n-1) subdiagonal elements of A.
   /// @param D Real array, size n.  On entry, D must contain the diagonal elements of A.
   /// @param DU Real array, size (n-1).  On entry, DU must contain the (n-1) superdiagonal elements of A.
   /// @ingroup NORM
   
   template <typename real_t>
   real_t LANGT(const char normType, const int_t n, complex<real_t> * const DL, complex<real_t> * const D, complex<real_t> * const DU)
   {
      const real_t zero = 0.0;
      const real_t one = 1.0;
      using std::isnan;
      using std::abs;
      real_t value = zero;
      if (n <= 0)
         return value;
      if (normType == 'M' || normType == 'm')
      {
         value = abs(D[n-1]);
         real_t temp;
         for (int_t i = 0; i < n-1; ++i)
         {
            temp = abs(DL[i]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
            temp = abs(D[i]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
            temp = abs(DU[i]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
         }
      }
      if (normType == 'O' || normType == 'o' || normType == '1')
      {
         if (n == 1)
            value = abs(D[0]);
         else
         {
            real_t temp;
            value = abs(D[0]) + abs(DL[0]);
            temp = abs(D[n-1]) + abs(DU[n-2]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
            for (int_t i = 1; i < n-1; ++i)
            {
               temp = abs(D[i])+abs(DU[i-1])+abs(DL[i]);
               if (value < temp)
                  value = temp;
               else if (isnan(temp))
                  return temp;
            }
         }
      }
      if (normType == 'I' || normType == 'i')
      {
         if (n == 1)
            value = abs(D[0]);
         else
         {
            real_t temp;
            value = abs(D[0]) + abs(DU[0]);
            temp = abs(D[n-1]) + abs(DL[n-2]);
            if (value < temp)
               value = temp;
            else if (isnan(temp))
               return temp;
            for (int_t i = 1; i < n-1; ++i)
            {
               temp = abs(D[i]) + abs(DL[i-1]) + abs(DU[i]);
               if (value < temp)
                  value = temp;
               else if (isnan(temp))
                  return temp;
            }
         }
      }
      if (normType == 'F' || normType == 'f' || normType == 'E' || normType == 'e')
      {
         real_t scale = zero;
         real_t sum = one;
      
         LATL::LASSQ(n, D, 1, scale, sum);
         if (n > 1)
         {
            LATL::LASSQ(n-1, DL, 1, scale, sum);
            LATL::LASSQ(n-1, DU, 1, scale, sum);
         }
         value = scale*sqrt(sum);
      }
      return value;
   }

}

#endif
