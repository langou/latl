//
//  lagtf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/31/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lagtf_h
#define _lagtf_h

/// @file lagtf.h Factorizes the matrix (T - lamda*I) where T is tridiagonal and lambda is a scalar.

#include <limits>
#include "latl.h"

namespace latl
{
   /// @brief Factorizes the matrix (T - lambda*I), where T is an n-by-n tridiagonal matrix and lambda is a scalar, as
   ///
   ///     T - lambda*I = P L U
   ///
   /// where P is a permutation matrix, L is a unit lower tridiagonal matrix with at most one non-zero subdiagonal elements per column and U is an upper triangular matrix with at most two non-zero super-diagonal elements per column.
   ///
   /// The factorization is obtained by Gaussian elimination with partial pivoting and implicit row scaling.
   ///
   /// The parameter lambda is included so that lagtf may be used in conjunction with lagts to obtain eigenvectors of T by inverse iteration.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix T.  n >= 0
   /// @param A Real array size n.  On entry, A should contain the diagonal elements of T.  On exit, A is overwritten by the n diagonal elements of the upper triangular matrix U.
   /// @param lambda Real scalar.
   /// @param B Real array size n-1.  On entry, B should contain the n-1 super-diagonal elements of T.  On exit, B is overwritten by the n-1 super-diagonal elements of the matrix U.
   /// @param C Real array size n-1.  On entry, C should contain the n-1 sub-diagonal elements of T.  On exit, C is overwritten by the n-1 sub-diagonal elements of the matrix L.
   /// @param tol Real scalar.  On entry, a relative tolerance used to indicate whether the matrix (T - lambda*I) is nearly singular.  Tol should be set to an approximation of the largest relative error in the elements of T.  For example, if the elements of T are correct to about 4 significant figures, tol should be set to 5*10^-4.  If the supplied tol is less than relative machine precision, then the machine eps will be used in place of tol.
   /// @param D Real array size n-2.  On exit, D is overwritten by the n-2 second super-diagonal elements of the matrix U.
   /// @param IN Integer array size n.  On exit, IN contains details of the permutation matrix P.  If an interchange occured at the kth step of the elimination, then IN[k] = 1; otherwise IN[k] = 0.  The final element of IN, IN[n-1], contains the smallest positive integer j such that the j-1th diagonal entry is less than or equal to the norm of T-lambda*I times tol.  If no such j exists, then IN[n-1] is returned as 0.  If such a j exists, then a diagonal element of U is small, indicating that T-lambda*I is singular or nearly singular.
   /// @ingroup TRF
   
   template< typename real_t >
   int lagtf(const int_t n, real_t * const A, const real_t lambda, real_t * const B, real_t * const C, real_t tol, real_t * const D, int_t * const IN)
   {
      using std::numeric_limits;

      if (n < 0)
         return -1;
      
      if (n == 0)
         return 0;
      
      A[0] -= lambda;
      IN[n-1] = 0;
      if (n == 1)
      {
         if (A[0] == 0)
            IN[n-1] = 1;
         return 0;
      }
      
      real_t const zero(0.0);
      real_t const eps = numeric_limits<real_t>::epsilon();
      tol = std::max(tol, eps);
      real_t scale1 = std::abs(A[0]) + std::abs(B[0]);
      real_t scale2, piv1, piv2, mult, temp;
      for (int_t k = 0; k < n-1; ++k)
      {
         A[k+1] -= lambda;
         scale2 = std::abs(C[k]) + std::abs(A[k+1]);
         if (k < n-2)
         {
            scale2 += std::abs(B[k+1]);
         }
         if (A[k] == 0)
            piv1 = zero;
         else
         {
            piv1 = std::abs(A[k]/scale1);
         }
         if (C[k] == 0)
         {
            IN[k] = 0;
            piv2 = 0;
            scale1 = scale2;
            if (k < n-2)
               D[k] = zero;
         }
         else
         {
            piv2 = std::abs(C[k])/scale2;
            if (piv2 <= piv1)
            {
               IN[k] = 0;
               scale1 = scale2;
               C[k] = C[k] /A[k];
               A[k+1] -= C[k]*B[k];
               if (k < n-2)
                  D[k] = zero;
            }
            else
            {
               IN[k] = 1;
               mult = A[k]/C[k];
               A[k] = C[k];
               temp = A[k+1];
               A[k+1] = B[k] - mult*temp;
               if (k < n-2)
               {
                  D[k] = B[k+1];
                  B[k+1] = -mult*D[k];
               }
               B[k] = temp;
               C[k] = mult;
            }
         }
         if (std::max(piv1, piv2) <= tol && IN[n-1] == 0)
            IN[n-1] = k+1;
      }
      if (std::abs(A[n-1]) <= scale1*tol && IN[n-1] == 0)
         IN[n-1] = n;
      return 0;
   }
}

#endif
