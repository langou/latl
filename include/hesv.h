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

namespace LATL
{
   /// @brief Computes the solution to a complex system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n Hermitian matrix and X and B are n-by-nrhs matrices.
   ///
   /// The Bunch-Kaufman diagonal pivoting method is used to factor A as
   ///
   ///         L^H*D*L, if uplo == 'L', or
   ///         U*D*U^H if uplo == 'U',
   ///
   /// where U (or L) is a product of permutation and unit upper/lower triangular matrices and D is Hermitian and block diagonal.
   ///
   /// The factored form of A is then used to solve the system of equations.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith pivot is exactly zero.  The factorization has been completed, but the block diagonal matrix D is
   /// exactly singular and division by zero will occur if it is used to solve a system of equations.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.
   /// @param A Complex array size ldA-by-n.  On entry, the Hermitian matrix A.  On exit, the block diagonal matrix D and the multipliers used to obtain the factor U or L.
   /// @param ldA Column length of the array A.
   /// @param ipiv Integer array size n.  On exit, contains the details of the interchanges of D.
   /// @param bsdv Bool array size n.  On exit, contains the details of the block structure of D.
   /// @param B Complex array size n-by-nrhs.  On exit, contains the solution X.
   /// @param ldB Column length of the matrix B.  ldB >= n
   /// @ingroup DRIV
   
   template<typename real_t>
   int_t HESV(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, int_t ldA, int_t * ipiv, bool * bsdv, complex<real_t> * const B, int_t ldB)
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
      
      int_t info = LATL::HETRF(uplo, n, A, ldA, ipiv, bsdv);
      if (info == 0)
      {
         info = LATL::HETRS(uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
      }
      return info;
   }

   /// @brief Computes the solution to a complex system of linear equations
   ///
   ///      A * X = B
   ///
   /// where A is an n-by-n Hermitian matrix and X and B are n-by-nrhs matrices.
   ///
   /// The Bunch-Kaufman diagonal pivoting method is used to factor A as
   ///
   ///         L^H*D*L, if uplo == 'L', or
   ///         U*D*U^H if uplo == 'U',
   ///
   /// where U (or L) is a product of permutation and unit upper/lower triangular matrices and D is Hermitian and block diagonal.
   ///
   /// The factored form of A is then used to solve the system of equations.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith pivot is exactly zero.  The factorization has been completed, but the block diagonal matrix D is
   /// exactly singular and division by zero will occur if it is used to solve a system of equations.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the symmetric matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.
   /// @param A Complex array size ldA-by-n.  On entry, the Hermitian matrix A.  On exit, the block diagonal matrix D and the multipliers used to obtain the factor U or L.
   /// @param ldA Column length of the array A.
   /// @param ipiv Integer array size n.  On exit, contains the details of the interchanges of D.
   /// @param bsdv Bool array size n.  On exit, contains the details of the block structure of D.
   /// @param B Complex array size n-by-nrhs.  On exit, contains the solution X.
   /// @param ldB Column length of the matrix B.  ldB >= n
   /// @param nb Block size for computing the factorization.
   /// @ingroup DRIV
   
   template<typename real_t>
   int_t HESV(const char uplo, const int_t n, const int_t nrhs, complex<real_t> * const A, int_t ldA, int_t * ipiv, bool * bsdv, complex<real_t> * const B, int_t ldB, int_t nb)
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
      
      int_t info = LATL::HETRF(uplo, n, A, ldA, ipiv, bsdv, nb);
      if (info == 0)
      {
         info = LATL::HETRS(uplo, n, nrhs, A, ldA, ipiv, bsdv, B, ldB);
      }
      return info;
   }
}

#endif
