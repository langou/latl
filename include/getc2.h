//
//  getc2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 6/21/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _getc2_h
#define _getc2_h

/// @file getc2.h Computes an LU factorization with complete pivoting of an n-by-n matrix.

#include <limits>
#include <cmath>
#include <algorithm>
#include "labad.h"
#include "swap.h"
#include "ger.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes an LU factorization with complete pivoting of an n-by-n real matrix.
   ///
   /// The factorization has the form
   /// A = P * L * U * Q
   /// where P and Q are permutation matrices.  L is lower triangular with unit diagonal
   /// elements and U is upper triangular.
   /// @return 0 if success.
   /// @return 1 if diagonal element was reset to minimum value to avoid singularity.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Real n-by-n matrix.  On entry, the n-by-n matrix to be factored.
   /// On exit, contains the factors L and U from the factorization.  The unit diagonal elements
   /// are not stored.
   /// @param ldA The column length of the array A.  ldA >= n
   /// @param ipiv Integer array size n.  Contains pivot indices; for 0 <= i <= n-1,
   /// row i of the matrix has been exchanged with row ipiv[i].
   /// @param jpiv Integer array size n.  Contains pivot indices; for 0 <= j <= n-1,
   /// column j of the matrix has been exchanged with column jpiv[j].
   
   template<typename real_t>
   int getc2(int_t n, real_t *A, int_t ldA, int_t *ipiv, int_t *jpiv)
   {
      using std::numeric_limits;
      using std::abs;
      using std::max;
      const real_t one(1.0);
      const real_t zero(0.0);
      const real_t eps = numeric_limits<real_t>::epsilon();

      real_t smlnum = numeric_limits<real_t>::min()/eps;
      real_t bignum = one/smlnum;
      labad(smlnum, bignum);

      int info=0;
      int_t ipv=0;
      int_t jpv=0;
      real_t smin;


      for(int_t i=0; i<n-1;i++)
      {
         real_t xmax=zero;
         for(int_t jp=i;jp<n;jp++)
         {
            for(int_t ip=i;ip<n;ip++)
            {
               real_t a=A[ip+jp*ldA];
               if(abs(a)>=xmax)
               {
                  xmax=abs(a);
                  ipv=ip;
                  jpv=jp;
               }
            }
         }
         if(i==0)
            smin=max(eps*xmax,smlnum);
         if(ipv!=i)
            swap(n,A+ipv,ldA,A+i,ldA);
         ipiv[i]=ipv;
         if(jpv!=i)
            swap(n,A+ldA*jpv,1,A+ldA*i,1);
         jpiv[i]=jpv;
         if(abs(A[i+i*ldA])<smin)
         {
            info=1;
            A[i+i*ldA]=smin;
         }
         for(int_t j=i+1;j<n;j++)
            A[j+i*ldA]=A[j+i*ldA]/A[i+i*ldA];
         ger(n-i-1,n-i-1,-one,A+i+1+i*ldA,1,A+i+(i+1)*ldA,ldA,A+i+1+(i+1)*ldA,ldA);
      }
      if(abs(A[n-1+(n-1)*ldA])<smin)
      {
         info=1;
         A[n-1+(n-1)*ldA]=smin;
      }
      jpiv[n-1]=ipiv[n-1]=n-1;
      return info;
   }
   
   /// @brief Computes an LU factorization with complete pivoting of an n-by-n complex matrix.
   ///
   /// The factorization has the form
   /// A = P * L * U * Q
   /// where P and Q are permutation matrices.  L is lower triangular with unit diagonal
   /// elements and U is upper triangular.
   /// @return 0 if success.
   /// @return 1 if diagonal element was reset to minimum value to avoid singularity.
   /// @tparam real_t Floating point type.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Complex n-by-n matrix.  On entry, the n-by-n matrix to be factored.
   /// On exit, contains the factors L and U from the factorization.  The unit diagonal elements
   /// are not stored.
   /// @param ldA The column length of the array A.  ldA >= n
   /// @param ipiv Integer array size n.  Contains pivot indices; for 0 <= i <= n-1,
   /// row i of the matrix has been exchanged with row ipiv[i].
   /// @param jpiv Integer array size n.  Contains pivot indices; for 0 <= j <= n-1,
   /// column j of the matrix has been exchanged with column jpiv[j].

   template<typename real_t>
   int getc2(int_t n, complex<real_t> *A, int_t ldA, int_t *ipiv, int_t *jpiv)
   {
      using std::numeric_limits;
      using std::abs;
      using std::max;
      const real_t one(1.0);
      const complex<real_t> minus_one(-1.0);
      const real_t zero(0.0);
      const real_t eps = numeric_limits<real_t>::epsilon();

      real_t smlnum = numeric_limits<real_t>::min()/eps;
      real_t bignum = one/smlnum;
      labad(smlnum, bignum);

      int info=0;
      int_t ipv=0;
      int_t jpv=0;
      real_t smin;


      for(int_t i=0; i<n-1;i++)
      {
         real_t xmax=zero;
         for(int_t jp=i;jp<n;jp++)
         {
            for(int_t ip=i;ip<n;ip++)
            {
               complex<real_t> a=A[ip+jp*ldA];
               if(abs(a)>=xmax)
               {
                  xmax=abs(a);
                  ipv=ip;
                  jpv=jp;
               }
            }
         }
         if(i==0)
            smin=max(eps*xmax,smlnum);
         if(ipv!=i)
            swap(n,A+ipv,ldA,A+i,ldA);
         ipiv[i]=ipv;
         if(jpv!=i)
            swap(n,A+ldA*jpv,1,A+ldA*i,1);
         jpiv[i]=jpv;
         if(abs(A[i+i*ldA])<smin)
         {
            info=1;
            A[i+i*ldA]=smin;
         }
         for(int_t j=i+1;j<n;j++)
            A[j+i*ldA]=A[j+i*ldA]/A[i+i*ldA];
         ger(n-i-1,n-i-1,minus_one,A+i+1+i*ldA,1,A+i+(i+1)*ldA,ldA,A+i+1+(i+1)*ldA,ldA);
      }
      if(abs(A[n-1+(n-1)*ldA])<smin)
      {
         info=1;
         A[n-1+(n-1)*ldA]=smin;
      }
      jpiv[n-1]=ipiv[n-1]=n-1;
      return info;
   }
}

#endif
