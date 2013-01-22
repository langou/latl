//
//  hetrf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/31/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _hetrf_h
#define _hetrf_h

/// @file hetrf.h Computes the factorization of a complex Hermitian matrix.

#include "hetf2.h"
#include "lahef.h"
#include "latl.h"

namespace latl
{
   /// @brief Computes the factorization of a complex Hermitian matrix A using the Bunch-Kaufman diagonal pivoting method.  The form of the factorization is
   ///
   ///         A = U*D*U^H if uplo = 'U'
   ///         A = L^H*D*L if uplo = 'L'
   ///
   /// where U (or L) is a product of permutation and unit upper (or lower) triangular matrices and D is Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
   ///
   /// This is the blocked version of the algorithm.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith pivot is exactly zero.  The factorization has been completed, but the block diagonal matrix D is exactly singular and division by zero will occur if it is used to solve a system of equations.
   /// @tparam real_t Floating point type.
   /// @param uplo Indicates whether the Hermitian matrix A is stored as upper triangular or lower triangular.  The other triangular part of A is not referenced.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Complex array size ldA-by-n.  On entry, the Hermitian matrix A.  On exit, the block diagonal matrix D and the multipliers used to obtain the factor U or L.
   /// @param ldA Column length of the array A.
   /// @param IPIV Integer array size n.  On exit, contains the details of the interchanges of D.
   /// @param BSDV Bool array size n.  On exit, contains the details of the block structure of D.  If BSDV[k] = 0, then rows and columns k and IPIV[k] were interchanged and D[k, k] is a 1-by-1 diagonal block.  If BSDV[k] = 1, then k is part of a 2-by-2 diagonal block.  In a 2 by 2 block, if uplo = 'U', and IPIV[k] = IPIV[k-1], then rows and columns k-1 and IPIV[k] were interchanged.  If uplo = 'L' and IPIV[k] = IPIV[k+1], then rows and columns k+1 and IPIV[k] were interchanged.
   /// @param nb Block size, optional.  Default value is 64.
   /// @ingroup TRF

   template< typename real_t >
   int_t hetrf(const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA, int_t * const IPIV, bool * const BSDV, const int_t nb = 64)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      int_t info = 0, k, temp = 0, kb;
      if (uplo == 'U' || uplo == 'u')
      {
         // in this section, k is not an index variable
         k = n;
         while ( k >= 1)
         {
            if ( k > nb)
            {
               temp = latl::lahef(uplo, k, nb, kb, A, ldA, IPIV, BSDV);
            }
            else
            {
               temp = latl::hetf2(uplo, k, A, ldA, IPIV, BSDV);
               kb = k;
            }
            
            if (info == 0 && temp != 0)
            {
               info = temp;
            }
            
            k -= kb;
         }
      }
      else
      {
         k = 0;
         complex<real_t> *Akk;
         while (k < n)
         {
            Akk = A+ldA*k+k;
            if ( k < n-nb)
            {
               temp = latl::lahef(uplo, n-k, nb, kb, Akk, ldA, IPIV+k, BSDV+k);
            }
            else
            {
               kb = n-k;
               temp = latl::hetf2(uplo, kb, Akk, ldA, IPIV+k, BSDV+k);
            }
            
            if (info == 0 && temp != 0)
            {
               info = temp+k;
            }
            
            for (int_t j = k; j < k+kb; ++j)
            {
               IPIV[j] += k;
            }
            
            k += kb;
         }
      }
      return info;
   }
}


#endif
