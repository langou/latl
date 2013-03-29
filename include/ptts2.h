//
//  ptts2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/20/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _ptts2_h
#define _ptts2_h

///@file ptts2.h Solves a system of linear equations A * X = B.

#include "scal.h"
#include "latl.h"

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
   /// @ingroup COMP

   template<typename real_t>
   int_t PTTS2(const int_t n, const int_t nrhs, real_t * const D, real_t * const E, real_t * const B, const int_t ldB)
   {
      if (n < 0)
         return -1;
      if (nrhs < 0)
         return -2;
      if (ldB < n)
         return -6;
      if (n == 0 || nrhs == 0)
         return 0;
      if (n == 1)
      {
         LATL::SCAL(nrhs, 1.0/D[0], B, ldB);
         return 0;
      }
      real_t * Bj = B;
      for (int_t j = 0; j < nrhs; ++j)
      {
         for (int_t i = 1; i < n; ++i)
         {
            Bj[i] -= Bj[i-1]*E[i-1];
         }
         Bj[n-1] /= D[n-1];
         for (int_t i = n-2; i >= 0; --i)
         {
            Bj[i] = Bj[i]/D[i] - Bj[i+1]*E[i];
         }
         Bj += ldB;
      }
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
   /// @param iuplo Indicated if the subdiagonal E should instead be processed as a superdiagonal.
   /// @param n Order of the matrix A.   n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param D Real array size n.  On entry, the n diagonal elements of the diagonal matrix D from the L*D*L^H factorization of A as computed by pttrf.
   /// @param E Complex array size (n-1).  On entry, the (n-1) subdiagonal elements of the unit bidiagonal factor L from the L*D*L^H factorization as computed by pttrf.
   /// @param B Complex array size ldB-by-nrhs.  On entry, the right hand side vectors B for the system of linear equations.  On exit, the solution vectors X.
   /// @param ldB Column length of the array B.
   /// @ingroup COMP

   template<typename real_t>
   int_t PTTS2(bool iuplo, const int_t n, const int_t nrhs, real_t * const D, complex<real_t> * const E, complex<real_t> * const B, const int_t ldB)
   {
      if ( n < 0)
         return -1;
      if (nrhs < 0)
         return -2;
      if (ldB < n)
         return -6;
      if (n == 0 || nrhs == 0)
         return 0;
      if (n == 1)
      {
         LATL::SCAL(nrhs, 1.0/D[0], B, ldB);
         return 0;
      }
      if (iuplo == 1)
      {
         complex<real_t> * Bj = B;
         for (int_t j = 0; j < nrhs; ++j)
         {
            for (int_t i = 1; i < n; ++i)
            {
               Bj[i] -= Bj[i-1]*conj(E[i-1]);
            }
            Bj[n-1] = Bj[n-1]/D[n-1];
            for (int_t i = n-2; i >= 0; --i)
            {
               Bj[i] = Bj[i]/D[i]-Bj[i+1]*E[i];
            }
            Bj += ldB;
         }
      }
      else
      {
         complex<real_t> * Bj = B;
         for (int_t j = 0; j < nrhs; ++j)
         {
            for (int_t i = 1; i < n; ++i)
            {
               Bj[i] -= Bj[i-1]*E[i-1];
            }
            Bj[n-1] = Bj[n-1]/D[n-1];
            for (int_t i = n-2; i >= 0; --i)
            {
               Bj[i] = Bj[i]/D[i]-Bj[i+1]*conj(E[i]);
            }
            Bj += ldB;
         }
      }
   }
}

#endif
