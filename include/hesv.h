//
//  hesv.h
//  Linear Algebra Template Library  
//
//  Created by Stephanie Patterson on 3/12/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _hesv_h
#define _hesv_h

/// @file hesv.h


#include "hetrf.h"
#include "hetrs.h"

namespace latl
{
   /// @brief Computes the solution to a complex system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n Hermitian matrix and X and B are n-by-nrhs matrices.
   ///
   /// The Bunch-Kaufman diagonal pivoting method is used to factor A as L^H*D*L, if uplo == 'L', or U*D*U^H if uplo == 'U',
   /// where U (or L) is a product of permutation and unit upper/lower triangular matrices and D is Hermitian and block diagonal.
   ///
   /// The factored form of A is then used to solve the system of equations.
   ///
   ///
   
   template<typename real_t>
   int_t hesv(const char uplo, const int_t n, const int_t nrhs, real_t * const A, int_t ldA, int_t * ipiv, bool * bsdv, real_t * const B, int_t ldB)
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
      
      int_t info = latl::hetrf(uplo, n, A, ldA, ipiv, bsdv);
      if (info == 0)
      {
         info = latl::hetrs(uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
      }
      return info;
   }
   
   template<typename real_t>
   int_t hesv(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, int_t ldA, int_t * ipiv, bool * bsdv, complex<real_t> * const B, int_t ldB)
   {
      return latl::hesv< complex<real_t> > (uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
   }
   
   template<typename real_t>
   int_t hesv(const char uplo, const int_t n, const int_t nrhs, real_t * const A, int_t ldA, int_t * ipiv, bool * bsdv, real_t * const B, int_t ldB, int_t nb = 32)
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
      
      int_t info = latl::hetrf(uplo, n, A, ldA, ipiv, bsdv, nb);
      if (info == 0)
      {
         info = latl::hetrs(uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
      }
      return info;
   }
   
   template<typename real_t>
   int_t hesv(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, int_t ldA, int_t * ipiv, bool * bsdv, complex<real_t> * const B, int_t ldB, int_t nb = 32)
   {
      return latl::hesv< complex<real_t> > (uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB, nb);
   }
}

#endif
