//
//  lantb.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 1/21/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lantb_h
#define _lantb_h

/// @file lantb.h Returns the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a triangular band matrix A.

#include "latl.h"
#include <cmath>
#include "lassq.h"

namespace LATL
{
   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a real triangular band matrix A.
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
   /// @param diag specifies whether or not A is unit triangular.  If diag = 'U' or 'u' then A is unit triangular; otherwise A is assumed to not be unit triangular.
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param k Specifies the number of super/sub-diagonals of the matrix A.  k >= 0
   /// @param AB Pointer to banded triangular n-by-n matrix A.  The bands of A are stored as rows,
   /// while preserving the columns of A.  If uplo = 'U' or 'u' then only the upper triangular part of A
   /// is referenced and the lower part is not referenced; the diagonal is stored on row k, the first super-diagonal
   /// in row k-1 starting in column 1, the second super-diagonal in row k-2 starting in column 2, and so on.
   /// As an example, consider the following upper triangular matrix with n=5 and k=2.  On the left is the matrix in
   /// standard upper triangular form, and on the right is the same matrix in upper triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h d a . . ]        [ . . a b c ]
   ///        [ . i e b . ]        [ . d e f g ]
   ///        [ . . j f c ]        [ h i j k l ]
   ///        [ . . . k g ]
   ///        [ . . . . l ]
   /// Otherwise only the lower triangular part of A is referenced and the upper part is not referenced;
   /// the diagonal is stored in row zero, the first sub-diagonal in row 1, and second sub-diagonal in row 2, and so on.
   /// As an example, consider the following lower triangular matrix with n=5 and k=2.  On the left is the matrix in standard
   /// lower triangular form, and on the right in the same matrix in lower triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h . . . . ]        [ h i j k l ]
   ///        [ d i . . . ]        [ d e f g . ]
   ///        [ a e j . . ]        [ a b c . . ]
   ///        [ . b f k . ]
   ///        [ . . c g l ]
   /// @param ldAB Column length of the matrix A.  ldAB >= m
   /// @param Work Workspace vector of length m (optional).  If not used, workspace will be allocated and
   /// deallocated internally; only used for the infinity norm.
   /// @ingroup NORM
   
   template <typename real_t>
   real_t LANTB(const char normType, const char uplo, const char diag, const int_t n, const int_t k, real_t * const AB, const int_t ldAB, real_t *Work=NULL)
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
         real_t temp(0.0);
         if (diag == 'U' || diag == 'u')
         {
            value = 1.0;
            if (uplo == 'U' || uplo == 'u')
            {
               ABj += ldAB;
               for (int_t j = 1; j < n; ++j)
               {
                  for (int_t i = max(k-j, intzero); i < k; ++i)
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
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 1; i < min(k+1, n-j); ++i)
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
         }
         else 
         {
            if (uplo == 'U' || uplo == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = max(k-j, intzero); i < k+1; ++i)
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
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < min(n-j, k+1); ++i)
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
         }
      }
      else if (normType == 'O' || normType == 'o' || normType == '1')
      {
         real_t * ABj = AB;
         real_t sum(0.0);
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (diag == 'U' || diag == 'u')
               {
                  sum = 1;
                  for (int_t i = max(k-j, intzero); i < k; ++i)
                  {
                     sum += abs(ABj[i]);
                  }
               }
               else
               {
                  sum = zero;
                  for (int_t i = max(k-j, intzero); i < k+1; ++i)
                  {
                     sum += abs(ABj[i]);
                  }
               }
               if (value < sum)
               {
                  value = sum;
               }
               else if (isnan(sum))
                  return sum;
               ABj+= ldAB;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (diag == 'U' || diag == 'u')
               {
                  sum = 1;
                  for (int_t i = 1; i < min(k+1, n-j); ++i)
                  {
                     sum += abs(ABj[i]);
                  }
               }
               else
               {
                  sum = zero;
                  for (int_t i = 0; i < min(k+1, n-j); ++i)
                  {
                     sum += abs(ABj[i]);
                  }
               }
               if (value < sum)
               {
                  value = sum;
               }
               else if (isnan(sum))
                  return sum;
               ABj+= ldAB;
            }
         }
      }
      else if ( normType == 'I' || normType == 'i')
      {
         bool allocate=(Work==NULL)?1:0;
         if(allocate)
            Work = new real_t[n];
         real_t * ABj = AB;
         real_t sum(0.0);
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Work[j] = 1.0;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = max(k-j, intzero); i < k; ++i)
                  {
                     Work[j+i-k] += abs(ABj[i]);
                  }
                  ABj += ldAB;
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Work[j] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = max(k-j, intzero); i < k+1; ++i)
                  {
                     Work[j+i-k] += abs(ABj[i]);
                  }
                  ABj += ldAB;
               }
            }
         }
         else 
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = 1.0;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 1; i < min(n-j, k+1); ++i)
                  {
                     Work[i+j] += abs(ABj[i]);
                  }
                  ABj += ldAB;
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
                  for (int_t i = 0; i < min(n-j, k+1); ++i)
                  {
                     Work[i+j] += abs(ABj[i]);
                  }
                  ABj += ldAB;
               }
            }
         }
         for (int_t i = 0; i < n; ++i)
         {
            sum = Work[i];
            if (sum > value)
            {
               value = sum;
            }
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
      else if ( normType == 'F' || normType == 'f' || normType == 'E' ||normType == 'e')
      {
         real_t scale(0.0), sum(1.0), * ABj = AB;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               scale = 1.0;
               sum = n;
               if (k > 0)
               {
                  ABj += ldAB+k-1;
                  for (int_t j = 1; j < n; ++j)
                  {
                     LATL::LASSQ(min(j, k), ABj, 1, scale, sum);
                     if (k-j > 0)
                        ABj += ldAB-1;
                     else
                        ABj += ldAB;
                  }
               }
            }
            else
            {
               scale = zero;
               sum = 1.0;
               ABj += k;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::LASSQ(min(j, k)+1, ABj, 1, scale, sum);
                  if (k-j > 0)
                     ABj += ldAB-1;
                  else
                     ABj += ldAB;
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               scale = 1.0;
               sum = n;
               ABj += 1;
               if (k > 0)
               {
                  for (int_t j = 0; j < n-1; ++j)
                  {
                     LATL::LASSQ(min(k, n-j-1), ABj, 1, scale, sum);
                     ABj += ldAB;
                  }
               }
            }
            else
            {
               scale = zero;
               sum = 1.0;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::LASSQ(min(k+1, n-j), ABj, 1, scale, sum);
                  ABj += ldAB;
               }
            }
         }
         value = scale*sqrt(sum);
      }
      return value;
   }

   /// @brief Calculates the value of the one norm, Frobenius norm, infinity norm, or element of largest absolute value of a complex triangular band matrix A.
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
   /// @param diag Specifies whether or not A is unit triangular.  If diag = 'U' or 'u' then A is unit triangular; otherwise A is assumed to not be unit triangular.
   /// @param n Number of columns to be included in the norm. n >= 0
   /// @param k Specifies the number of super/sub-diagonals of the matrix A.  k >= 0
   /// @param AB Pointer to banded triangular n-by-n matrix A.  The bands of A are stored as rows,
   /// while preserving the columns of A.  If uplo = 'U' or 'u' then only the upper triangular part of A
   /// is referenced and the lower part is not referenced; the diagonal is stored on row k, the first super-diagonal
   /// in row k-1 starting in column 1, the second super-diagonal in row k-2 starting in column 2, and so on.
   /// As an example, consider the following upper triangular matrix with n=5 and k=2.  On the left is the matrix in
   /// standard upper triangular form, and on the right is the same matrix in upper triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h d a . . ]        [ . . a b c ]
   ///        [ . i e b . ]        [ . d e f g ]
   ///        [ . . j f c ]        [ h i j k l ]
   ///        [ . . . k g ]
   ///        [ . . . . l ]
   /// Otherwise only the lower triangular part of A is referenced and the upper part is not referenced;
   /// the diagonal is stored in row zero, the first sub-diagonal in row 1, and second sub-diagonal in row 2, and so on.
   /// As an example, consider the following lower triangular matrix with n=5 and k=2.  On the left is the matrix in standard
   /// lower triangular form, and on the right in the same matrix in lower triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h . . . . ]        [ h i j k l ]
   ///        [ d i . . . ]        [ d e f g . ]
   ///        [ a e j . . ]        [ a b c . . ]
   ///        [ . b f k . ]
   ///        [ . . c g l ]
   /// @param ldAB Column length of the matrix A.  ldAB >= m
   /// @param Work Workspace vector of length m (optional).  If not used, workspace will be allocated and
   /// deallocated internally; only used for the infinity norm.
   /// @ingroup NORM
   
   template <typename real_t>
   real_t LANTB(const char normType, const char uplo, const char diag, const int_t n, const int_t k, complex<real_t> * const AB, const int_t ldAB, real_t *Work=NULL)
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
         real_t temp(0.0);
         if (diag == 'U' || diag == 'u')
         {
            value = 1.0;
            if (uplo == 'U' || uplo == 'u')
            {
               ABj += ldAB;
               for (int_t j = 1; j < n; ++j)
               {
                  for (int_t i = max(k-j, intzero); i < k; ++i)
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
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 1; i < min(k+1, n-j); ++i)
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
         }
         else 
         {
            if (uplo == 'U' || uplo == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = max(k-j, intzero); i < k+1; ++i)
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
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 0; i < min(n-j, k+1); ++i)
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
         }
      }
      else if (normType == 'O' || normType == 'o' || normType == '1')
      {
         complex<real_t> * ABj = AB;
         real_t sum(0.0);
         if (uplo == 'U' || uplo == 'u')
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (diag == 'U' || diag == 'u')
               {
                  sum = 1;
                  for (int_t i = max(k-j, intzero); i < k; ++i)
                  {
                     sum += abs(ABj[i]);
                  }
               }
               else
               {
                  sum = zero;
                  for (int_t i = max(k-j, intzero); i < k+1; ++i)
                  {
                     sum += abs(ABj[i]);
                  }
               }
               if (value < sum)
               {
                  value = sum;
               }
               else if (isnan(sum))
                  return sum;
               ABj+= ldAB;
            }
         }
         else
         {
            for (int_t j = 0; j < n; ++j)
            {
               if (diag == 'U' || diag == 'u')
               {
                  sum = 1;
                  for (int_t i = 1; i < min(k+1, n-j); ++i)
                  {
                     sum += abs(ABj[i]);
                  }
               }
               else
               {
                  sum = zero;
                  for (int_t i = 0; i < min(k+1, n-j); ++i)
                  {
                     sum += abs(ABj[i]);
                  }
               }
               if (value < sum)
               {
                  value = sum;
               }
               else if (isnan(sum))
                  return sum;
               ABj+= ldAB;
            }
         }
      }
      else if ( normType == 'I' || normType == 'i')
      {
         bool allocate=(Work==NULL)?1:0;
         if(allocate)
            Work = new real_t[n];
         complex<real_t> * ABj = AB;
         real_t sum(0.0);
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Work[j] = 1.0;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = max(k-j, intzero); i < k; ++i)
                  {
                     Work[j+i-k] += abs(ABj[i]);
                  }
                  ABj += ldAB;
               }
            }
            else
            {
               for (int_t j = 0; j < n; ++j)
               {
                  Work[j] = zero;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = max(k-j, intzero); i < k+1; ++i)
                  {
                     Work[j+i-k] += abs(ABj[i]);
                  }
                  ABj += ldAB;
               }
            }
         }
         else 
         {
            if (diag == 'U' || diag == 'u')
            {
               for (int_t i = 0; i < n; ++i)
               {
                  Work[i] = 1.0;
               }
               for (int_t j = 0; j < n; ++j)
               {
                  for (int_t i = 1; i < min(n-j, k+1); ++i)
                  {
                     Work[i+j] += abs(ABj[i]);
                  }
                  ABj += ldAB;
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
                  for (int_t i = 0; i < min(n-j, k+1); ++i)
                  {
                     Work[i+j] += abs(ABj[i]);
                  }
                  ABj += ldAB;
               }
            }
         }
         for (int_t i = 0; i < n; ++i)
         {
            sum = Work[i];
            if (sum > value)
            {
               value = sum;
            }
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
      else if ( normType == 'F' || normType == 'f' || normType == 'E' ||normType == 'e')
      {
         real_t scale(0.0), sum(1.0);
         complex<real_t> * ABj = AB;
         if (uplo == 'U' || uplo == 'u')
         {
            if (diag == 'U' || diag == 'u')
            {
               scale = 1.0;
               sum = n;
               if (k > 0)
               {
                  ABj += ldAB+k-1;
                  for (int_t j = 1; j < n; ++j)
                  {
                     LATL::LASSQ(min(j, k), ABj, 1, scale, sum);
                     if (k-j > 0)
                        ABj += ldAB-1;
                     else
                        ABj += ldAB;
                  }
               }
            }
            else
            {
               scale = zero;
               sum = 1.0;
               ABj += k;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::LASSQ(min(j, k)+1, ABj, 1, scale, sum);
                  if (k-j > 0)
                     ABj += ldAB-1;
                  else
                     ABj += ldAB;
               }
            }
         }
         else
         {
            if (diag == 'U' || diag == 'u')
            {
               scale = 1.0;
               sum = n;
               ABj += 1;
               if (k > 0)
               {
                  for (int_t j = 0; j < n-1; ++j)
                  {
                     LATL::LASSQ(min(k, n-j-1), ABj, 1, scale, sum);
                     ABj += ldAB;
                  }
               }
            }
            else
            {
               scale = zero;
               sum = 1.0;
               for (int_t j = 0; j < n; ++j)
               {
                  LATL::LASSQ(min(k+1, n-j), ABj, 1, scale, sum);
                  ABj += ldAB;
               }
            }
         }
         value = scale*sqrt(sum);
      }
      return value;
   }

}

#endif
