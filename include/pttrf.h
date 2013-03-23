//
//  pttrf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/22/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _pttrf_h
#define _pttrf_h

/// @file pttrf.h Computes the L*D*L' factorization of a symmetric positive definite tridiagonal matrix A.

#include "latl.h"

namespace LATL
{
   /// @brief Computes the L*D*L' factorization of a real symmetric positive definite tridiagonal matrix A.  The factorization may also be regarded as having the form A = U'*D*U.
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
   int_t pttrf(const int_t n, real_t * const D, real_t * const E)
   {
      if (n < 0)
         return -1;
      
      if ( n == 0)
         return 0;
      
      const real_t zero(0.0);
      int_t i4;
      real_t ei;
      
      i4 = (n-1)%4;
      for (int_t i = 0; i < i4; ++i)
      {
         if (D[i] <= zero)
            return i+1;
         ei = E[i];
         E[i] = ei/D[i];
         D[i+1] -= E[i]*ei;
      }
      
      for (int_t i = i4; i < n-4; i += 4)
      {
         
         if (D[i] <= zero)
         {
            return i+1;
         }
         
         ei = E[i];
         E[i] = ei/D[i];
         D[i+1] -= E[i]*ei;
         
         if (D[i+1] <= zero)
            return i+2;
         
         ei = E[i+1];
         E[i+1] = ei/D[i+1];
         D[i+2] -= E[i+1] * ei;
         
         
         if (D[i+2] <= zero)
            return i+3;
         
         ei = E[i+2];
         E[i+2] = ei/D[i+2];
         D[i+3] -= E[i+2] * ei;
         
         if (D[i+3] <= zero)
            return i+4;
         
         ei = E[i+3];
         E[i+3] = ei/D[i+3];
         D[i+4] -= E[i+3] * ei;
      }
      
      if (D[n-1] <= zero)
         return n;
      
      return 0;
   }
   
   /// @brief Computes the L*D*L^H factorization of a complex Hermitian positive definite tridiagonal matrix A.  The factorization may also be regarded as having the form A = U^H*D*U.
   ///
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the leading minor of order i is not positive definite.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A.  n >= 0
   /// @param D Real array, size n.  On entry, the n diagonal elements of the tridiagonal matrix A.  On exit, the n diagonal elements of the diagonal matrix D from the L*D*L^H factorization of A.
   /// @param E Complex array, size (n-1).  On entry, the (n-1) subdiagonal elements of the tridiagonal matrix A.  On exit, the (n-1) subdiagonal elements of the unit bidiagonal factor L from the L*D*L' factorization of A.  E can also be regarded as the superdiagonal of the unit bidiagonal factor U from the U'*D*U factorization of A.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t pttrf(const int_t n, real_t * const D, complex<real_t> * const E)
   {
      if (n < 0)
         return -1;
      
      if ( n == 0)
         return 0;
      
      const real_t zero(0.0);
      int_t i4 = (n-1)%4;
      real_t eir, eii, f, g;
      
      for (int_t i = 0; i < i4; ++i)
      {
         if (D[i] <= zero)
            return i+1;
         eir = real(E[i]);
         eii = imag(E[i]);
         f = eir/D[i];
         g = eii/D[i];
         E[i] = complex<real_t>(f, g);
         D[i+1] -= f*eir + g*eii;
      }
      
      for (int_t i = i4; i < n-4; i += 4)
      {
         if (D[i] <= zero)
         {
            return i+1;
         }
         
         eir = real(E[i]);
         eii = imag(E[i]);
         f = eir/D[i];
         g = eii/D[i];
         E[i] = complex<real_t>(f, g);
         D[i+1] -= f*eir + g*eii;
         
         if (D[i+1] <= zero)
            return i+2;
         
         eir = real(E[i+1]);
         eii = imag(E[i+1]);
         f = eir/D[i+1];
         g = eii/D[i+1];
         E[i+1] = complex<real_t>(f, g);
         D[i+2] -= f*eir + g*eii;
         
         
         if (D[i+2] <= zero)
            return i+3;
         
         eir = real(E[i+2]);
         eii = imag(E[i+2]);
         f = eir/D[i+2];
         g = eii/D[i+2];
         E[i+2] = complex<real_t>(f, g);
         D[i+3] -= f*eir + g*eii;
         
         if (D[i+3] <= zero)
            return i+4;
         
         eir = real(E[i+3]);
         eii = imag(E[i+3]);
         f = eir/D[i+3];
         g = eii/D[i+3];
         E[i+3] = complex<real_t>(f, g);
         D[i+4] -= f*eir + g*eii;
      }
      
      if (D[n-1] <= zero)
         return n;
      
      return 0;
   }
}

#endif
