//
//  pttrs.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/21/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _pttrs_h
#define _pttrs_h

/// @file pttrs.h Solves a system of linear equations A * X = B.

#include "latl.h"
#include "ptts2.h"

namespace LATL
{
   /// @brief Solves a system of linear equations A * X = B.
   ///
   ///      A * X = B
   ///
   /// where A is a tridiagonal matrix with L*D*L' or U'*D*U factorization previously computed by pttrf.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A.   n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param D Real array size n.  On entry, the n diagonal elements of the diagonal matrix D from the L*D*L' factorization of A as computed by pttrf.
   /// @param E Real array size (n-1).  On entry, the (n-1) subdiagonal elements of the unit bidiagonal factor L from the L*D*L' factorization as computed by pttrf.
   /// @param B Real array size ldB-by-nrhs.  On entry, the right hand side vectors B for the system of linear equations.  On exit, the solution vectors X.
   /// @param ldB Column length of the array B.
   /// @param nb Block size.
   
   template< typename real_t>
   int_t PTTRS(const int_t n, const int_t nrhs, real_t * const D, real_t * const E, real_t * const B, const int_t ldB, const int_t nb)
   {
      if ( n < 0)
         return -1;
      if ( nrhs < 0)
         return -2;
      if (ldB < n)
         return -6;
      
      if ( n == 0 || nrhs == 0)
         return 0;
      using std::min;
      if (nb >= nrhs)
         LATL::PTTS2(n, nrhs, D, E, B, ldB);
      else
      {
         real_t * Bj = B;
         for (int_t j = 0; j < nrhs; j += nb)
         {
            int_t jb = min(nrhs-j, nb);
            LATL::PTTS2(n, jb, D, E, Bj, ldB);
            Bj += ldB*nb;
         }
      }
      return 0;
   }
   
   /// @brief Solves a system of linear equations A * X = B.
   ///
   ///      A * X = B
   ///
   /// where A is a Hermitian tridiagonal matrix with L*D*L^H or U^H*D*U factorization previously computed by pttrf.
   ///
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the array E is a subdiagonal or superdiagonal.
   /// @param n Order of the matrix A.   n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param D Real array size n.  On entry, the n diagonal elements of the diagonal matrix D from the L*D*L^H factorization of A as computed by pttrf.
   /// @param E Complex array size (n-1).  On entry, the (n-1) subdiagonal elements of the unit bidiagonal factor L from the L*D*L^H factorization as computed by pttrf.
   /// @param B Complex array size ldB-by-nrhs.  On entry, the right hand side vectors B for the system of linear equations.  On exit, the solution vectors X.
   /// @param ldB Column length of the array B.
   /// @param nb Block size.
   
   template< typename real_t>
   int_t PTTRS(const char uplo, const int_t n, const int_t nrhs, real_t * const D, complex<real_t> * const E, complex<real_t> * const B, const int_t ldB, const int_t nb)
   {
      if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldB < n)
         return -7;
      
      using std::min;
      if (n == 0 || nrhs == 0)
         return 0;
      bool iuplo = 0;
      if (uplo == 'U' || uplo == 'u')
         iuplo = 1;
      if (nb > nrhs)
         LATL::PTTS2(iuplo, n, nrhs, D, E, B, ldB);
      else
      {
         complex<real_t> * Bj = B;
         for (int_t j = 0; j < nrhs; j += nb)
         {
            int_t jb = min(nrhs-j, nb);
            LATL::PTTS2(iuplo, n, jb, D, E, Bj, ldB);
            Bj += ldB*nb;
         }
      }
      return 0;
   }
}

#endif
