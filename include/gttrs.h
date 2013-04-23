//
//  gttrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/22/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _gttrs_h
#define _gttrs_h

#include "latl.h"
#include "gtts2.h"

namespace LATL
{
   /// @brief Solves a real system of equations of the form
   ///
   ///      A * X = B      or    A' * X = B
   ///
   /// using the LU factorization of the tridiagonal matrix A computed by GTTRF.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies the form of the system of equations.
   ///
   ///      'N' = no transpose
   ///      'T' or 'C' = transpose
   ///
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param DL Real array, size n-1.  On entry, should contain the (n-1) multipliers that define the matrix L from the LU factorization of A.
   /// @param D Real array, size n.  On entry, should contain the n diagonal elements of the upper triangular matrix U from the LU factorization of A.
   /// @param DU Real array, size n-1.  On entry, the (n-1) elements of the first superdiagonal of U.
   /// @param DU2 Real array, size n-2.  On entry, the (n-2) elements of the second superdiagonal of U.
   /// @param ipiv Integer array, size n.  The pivot indices from the factorization LU.
   /// @param B Real matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B.
   /// @param nb Block size.
   /// @ingroup COMP
   
   template< typename real_t>
   int_t GTTRS(const char trans, const int_t n, const int_t nrhs, real_t * const DL, real_t * const D, real_t * const DU, real_t * const DU2, int_t * const ipiv, real_t * const B, const int_t ldB, const int_t nb)
   {
      if (trans != 'N' && trans != 'T' && trans != 'C' && trans != 'n' && trans != 't' && trans != 'c')
         return -1;
      if ( n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldB < n)
         return -10;
      if (n == 0 || nrhs == 0)
         return 0;
      int_t itrans = 1;
      using std::min;
      if (trans == 'N' || trans == 'n')
         itrans = 0;
      if (nb >= nrhs)
         LATL::GTTS2(itrans, n, nrhs, DL, D, DU, DU2, ipiv, B, ldB);
      else
      {
         real_t * Bj = B;
         for (int_t j = 0; j < nrhs; j += nb)
         {
            int_t jb = min(nrhs-j, nb);
            LATL::GTTS2(itrans, n, jb, DL, D, DU, DU2, ipiv, Bj, ldB);
            Bj += ldB*nb;
         }
      }
      return 0;
   }
   
   /// @brief Solves a complex system of equations of the form
   ///
   ///      A * X = B      or    A.' * X = B     or    A' * X = B
   ///
   /// using the LU factorization of the tridiagonal matrix A computed by GTTRF.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies the form of the system of equations.
   ///
   ///      'N' = no transpose
   ///      'T' = transpose
   ///      'C' = conjugate transpose
   ///
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param DL Complex array, size n-1.  On entry, should contain the (n-1) multipliers that define the matrix L from the LU factorization of A.
   /// @param D Complex array, size n.  On entry, should contain the n diagonal elements of the upper triangular matrix U from the LU factorization of A.
   /// @param DU Complex array, size n-1.  On entry, the (n-1) elements of the first superdiagonal of U.
   /// @param DU2 Complex array, size n-2.  On entry, the (n-2) elements of the second superdiagonal of U.
   /// @param ipiv Integer array, size n.  The pivot indices from the factorization LU.
   /// @param B Complex matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B.
   /// @param nb Block size.
   /// @ingroup COMP
   
   template< typename real_t>
   int_t GTTRS(const char trans, const int_t n, const int_t nrhs, complex<real_t> * const DL, complex<real_t> * const D, complex<real_t> * const DU, complex<real_t> * const DU2, int_t * const ipiv, complex<real_t> * const B, const int_t ldB, const int_t nb)
   {
      if (trans != 'N' && trans != 'T' && trans != 'C' && trans != 'n' && trans != 't' && trans != 'c')
         return -1;
      if ( n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldB < n)
         return -10;
      if (n == 0 || nrhs == 0)
         return 0;
      int_t itrans = 0;
      using std::min;
      if (trans == 'T' || trans == 't')
         itrans = 1;
      if (trans == 'C' || trans == 'c')
         itrans = 2;
      if (nb >= nrhs)
         LATL::GTTS2(itrans, n, nrhs, DL, D, DU, DU2, ipiv, B, ldB);
      else
      {
         complex<real_t> * Bj = B;
         for (int_t j = 0; j < nrhs; j += nb)
         {
            int_t jb = min(nrhs-j, nb);
            LATL::GTTS2(itrans, n, jb, DL, D, DU, DU2, ipiv, Bj, ldB);
            Bj += ldB*nb;
         }
      }
      return 0;
   }
}


#endif
