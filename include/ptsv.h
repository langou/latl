//
//  ptsv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/21/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _ptsv_h
#define _ptsv_h

/// @file ptsv.h Computes the solution to a system of linear equations A*X=B.

#include "latl.h"
#include "pttrs.h"

namespace LATL
{
   /// @brief Computes the solution to a real system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n tridiagonal symmetric matrix and X and B are n-by-nrhs matrices.
   ///
   /// The L*D*L' factorization of A, first computed by pttrf.h, is then used to solve the system of equations.
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.  
   /// @param n Order of the matrix A.  n >= 0
   /// @param D Real array, size n.  On entry, the n diagonal elements of the tridiagonal matrix A.  On exit, the n diagonal elements of the diagonal matrix D from the L*D*L' factorization of A.
   /// @param E Real array, size (n-1).  On entry, the (n-1) subdiagonal elements of the tridiagonal matrix A.  On exit, the (n-1) subdiagonal elements of the unit bidiagonal factor L from the L*D*L' factorization of A.  E can also be regarded as the superdiagonal of the unit bidiagonal factor U from the U'*D*U factorization of A.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t PTSV(const int_t n, const int_t nrhs, real_t * const D, real_t * const E, real_t * const B, const int_t ldB, const int_t nb)
   {
      if ( n < 0)
         return -1;
      if (nrhs < 0)
         return -2;
      if (ldB<n)
         return -6;
      
      int_t info;
      info = LATL::PTTRF(n, D, E);
      if (info == 0)
      {
         info = LATL::PTTRS(n, nrhs, D, E, B, ldB, nb);
      }
      return info;
   }
   
   /// @brief Computes the solution to a complex system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n tridiagonal Hermitian matrix and X and B are n-by-nrhs matrices.
   ///
   /// The L*D*L^H factorization of A, first computed by pttrf.h, is then used to solve the system of equations.
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A.  n >= 0
   /// @param D Real array, size n.  On entry, the n diagonal elements of the tridiagonal matrix A.  On exit, the n diagonal elements of the diagonal matrix D from the L*D*L^H factorization of A.
   /// @param E Complex array, size (n-1).  On entry, the (n-1) subdiagonal elements of the tridiagonal matrix A.  On exit, the (n-1) subdiagonal elements of the unit bidiagonal factor L from the L*D*L^H factorization of A.  E can also be regarded as the superdiagonal of the unit bidiagonal factor U from the U^H*D*U factorization of A.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t PTSV(const int_t n, const int_t nrhs, real_t * const D, complex<real_t> * E, complex<real_t> * B, const int_t ldB, const int_t nb)
   {
      if ( n < 0)
         return -1;
      if (nrhs < 0)
         return -2;
      if (ldB<n)
         return -6;
      
      int_t info;
      info = LATL::PTTRF(n, D, E);
      if (info == 0)
      {
         info = LATL::PTTRS('L', n, nrhs, D, E, B, ldB, nb);
      }
      return info;
   }
}


#endif
