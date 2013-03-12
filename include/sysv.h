//
//  sysv.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 2/19/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _sysv_h
#define _sysv_h

/// @file sysv.h Computes the solution to a system of linear equations A*X = B.

#include "sytrf.h"
#include "sytrs.h"

namespace latl
{
   /// @brief Computes the solution to a real system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n symmetric matrix and X and B are n-by-nrhs matrices.
   ///
   /// The Bunch-Kaufman diagonal pivoting method is used to factor A as L'*D*L, if uplo == 'L', or U*D*U' if uplo == 'U',
   /// where U (or L) is a product of permutation and unit upper/lower triangular matrices and D is symmetric and block diagonal.
   ///
   /// The factored form of A is then used to solve the system of equations.
   ///
   ///
   
   template<typename real_t>
   int_t sysv(const char uplo, const int_t n, const int_t nrhs, real_t * const A, int_t ldA, int_t * ipiv, bool * bsdv, real_t * const B, int_t ldB)
   {
      if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldA < n)
         return -5;
      if (ldB < n)
         return -8;
      
      int_t info = latl::sytrf(uplo, n, A, ldA, ipiv, bsdv);
      if (info == 0)
      {
         info = latl::sytrs(uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
      }
      return info;
   }
   
   template<typename real_t>
   int_t sysv(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, int_t ldA, int_t * ipiv, bool * bsdv, complex<real_t> * const B, int_t ldB)
   {
      return latl::sysv< complex<real_t> > (uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
   }
   
   template<typename real_t>
   int_t sysv(const char uplo, const int_t n, const int_t nrhs, real_t * const A, int_t ldA, int_t * ipiv, bool * bsdv, real_t * const B, int_t ldB, int_t nb = 32)
   {
      if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldA < n)
         return -5;
      if (ldB < n)
         return -8;
      
      int_t info = latl::sytrf(uplo, n, A, ldA, ipiv, bsdv, nb);
      if (info == 0)
      {
         info = latl::sytrs(uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
      }
      return info;
   }
   
   template<typename real_t>
   int_t sysv(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, int_t ldA, int_t * ipiv, bool * bsdv, complex<real_t> * const B, int_t ldB, int_t nb = 32)
   {
      return latl::sysv< complex<real_t> > (uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB, nb);
   }
}
#endif
