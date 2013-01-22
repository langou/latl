//
//  gttrf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/22/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _gttrf_h
#define _gttrf_h

/// @file gttrf.h Computes an LU factorization of a tridiagonal matrix A.

#include "latl.h"

namespace latl
{
   /// @brief Computes an LU factorization of a real tridiagonal matrix A using elimination with partial pivoting and row interchanges.
   ///
   /// The factorization has the form
   ///
   ///         A = L * U
   ///
   /// where L is a product of permutation and unit lower bidiagonal matrices and U is upper triangular with nonzeros only in the main diagonal and first two superdiagonals.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith diagonal element is exactly zero.  The factorization has been completed, but the factor U is exactly singular, and division by zero will occur if it is used to solve a system of equations.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A.  n >= 0
   /// @param DL Real array, size n-1.  On entry, DL must contain the (n-1) subdiagonal elements of A.  On exit, DL is overwritten by the (n-1) multipliers that define the matrix L from the LU factorization of A.
   /// @param D Real array, size n.  On entry, D must contain the diagonal elements of A.  On exit, D is overwritten by the n diagonal elements of the upper triangular matrix U from the LU factorization of A.
   /// @param DU Real array, size (n-1).  On entry, DU must contain the (n-1) superdiagonal elements of A.  On exit, DU is overwritten by the (n-1) elements of the first superdiagonal of U.
   /// @param DU2 Real array, size (n-2).  On exit, DU2 is overwritten by the n-2 elements of the second superdiagonal of U.
   /// @param IPIV Integer array, size n.  The pivot indices; for 0 <= i < n, row i of the matrix was interchanged with row IPIV[i].  IPIV[i] will always be either i or i+1; IPIV[i] = i indicates a row interchange was not required.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t gttrf(const int_t n, real_t * const DL, real_t * const D, real_t * const DU, real_t * const DU2, int_t * const IPIV)
   {
      if (n < 0)
         return -1;
      
      if (n == 0)
         return 0;
      
      const real_t zero(0.0);
      for (int_t i = 0; i < n; ++i)
         IPIV[i] = i;
      for (int_t i = 0; i < n-2; ++i)
         DU2[i] = zero;
      
      real_t temp;
      for (int_t i = 0; i < n-2; ++i)
      {
         if (std::abs(D[i]) >= std::abs(DL[i]))
         {
            if (D[i] != zero)
            {
               DL[i] = DL[i]/D[i];
               D[i+1] -= DL[i] *DU[i];
            }
         }
         else
         {
            temp = D[i]/DL[i];
            D[i] = DL[i];
            DL[i] = temp;
            temp = DU[i];
            DU[i] = D[i+1];
            D[i+1] = temp - DL[i]*D[i+1];
            DU2[i] = DU[i+1];
            DU[i+1] = -DL[i]*DU[i+1];
            ++IPIV[i];
         }
      }
      if (n > 1)
      {
         int_t i = n-2;
         if (std::abs(D[i]) >= std::abs(DL[i]))
         {
            if (D[i] != zero)
            {
               DL[i] = DL[i] / D[i];
               D[i+1] -= DL[i] *DU[i];
            }
         }
         else
         {
            temp = D[i]/DL[i];
            D[i] = DL[i];
            DL[i] = temp;
            temp = DU[i];
            DU[i] = D[i+1];
            D[i+1] = temp - DL[i]*D[i+1];
            ++IPIV[i];
         }
      }
      for (int_t i = 0; i < n; ++i)
      {
         if (D[i] == zero)
            return i+1;
      }
      return 0;
   }
   
   /// @brief Computes an LU factorization of a complex tridiagonal matrix A using elimination with partial pivoting and row interchanges.
   ///
   /// The factorization has the form
   ///
   ///         A = L * U
   ///
   /// where L is a product of permutation and unit lower bidiagonal matrices and U is upper triangular with nonzeros only in the main diagonal and first two superdiagonals.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith diagonal element is exactly zero.  The factorization has been completed, but the factor U is exactly singular, and division by zero will occur if it is used to solve a system of equations.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A.  n >= 0
   /// @param DL Complex array, size n-1.  On entry, DL must contain the (n-1) subdiagonal elements of A.  On exit, DL is overwritten by the (n-1) multipliers that define the matrix L from the LU factorization of A.
   /// @param D Complex array, size n.  On entry, D must contain the diagonal elements of A.  On exit, D is overwritten by the n diagonal elements of the upper triangular matrix U from the LU factorization of A.
   /// @param DU Complex array, size (n-1).  On entry, DU must contain the (n-1) superdiagonal elements of A.  On exit, DU is overwritten by the (n-1) elements of the first superdiagonal of U.
   /// @param DU2 Complex array, size (n-2).  On exit, DU2 is overwritten by the n-2 elements of the second superdiagonal of U.
   /// @param IPIV Integer array, size n.  The pivot indices; for 0 <= i < n, row i of the matrix was interchanged with row IPIV[i].  IPIV[i] will always be either i or i+1; IPIV[i] = i indicates a row interchange was not required.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t gttrf(const int_t n, complex<real_t> * const DL, complex<real_t> * const D, complex<real_t> * const DU, complex<real_t> * const DU2, int_t * const IPIV)
   {
      if (n < 0)
         return -1;
      
      if (n == 0)
         return 0;
      
      const real_t zero(0.0);
      for (int_t i = 0; i < n; ++i)
         IPIV[i] = i;
      for (int_t i = 0; i < n-2; ++i)
         DU2[i] = zero;
      
      complex<real_t> temp;
      for (int_t i = 0; i < n-2; ++i)
      {
         if (std::abs(D[i]) >= std::abs(DL[i]))
         {
            if (std::abs(D[i]) != zero)
            {
               DL[i] = DL[i]/D[i];
               D[i+1] -= DL[i] *DU[i];
            }
         }
         else
         {
            temp = D[i]/DL[i];
            D[i] = DL[i];
            DL[i] = temp;
            temp = DU[i];
            DU[i] = D[i+1];
            D[i+1] = temp - DL[i]*D[i+1];
            DU2[i] = DU[i+1];
            DU[i+1] = -DL[i]*DU[i+1];
            ++IPIV[i];
         }
      }
      if (n > 1)
      {
         int_t i = n-2;
         if (std::abs(D[i]) >= std::abs(DL[i]))
         {
            if (std::abs(D[i]) != zero)
            {
               DL[i] = DL[i] / D[i];
               D[i+1] -= DL[i] *DU[i];
            }
         }
         else
         {
            temp = D[i]/DL[i];
            D[i] = DL[i];
            DL[i] = temp;
            temp = DU[i];
            DU[i] = D[i+1];
            D[i+1] = temp - DL[i]*D[i+1];
            ++IPIV[i];
         }
      }
      for (int_t i = 0; i < n; ++i)
      {
         if (std::abs(D[i]) == zero)
            return i+1;
      }
      return 0;
   }
}

#endif
