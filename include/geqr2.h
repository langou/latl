//
//  geqr2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/26/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _geqr2_h
#define _geqr2_h

/// @file geqr2.h Computes a QR factorization of a matrix.

#include <algorithm>
#include "larf.h"
#include "larfg.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes a QR factorization of a real matrix A.
   ///
   /// The matrix Q is represented as a product of elementary reflectors
   ///
   ///         Q = H_1 H_2 . . . H_k, where k = min(m,n).
   /// Each H_i has the form
   ///
   ///         H_i = I - tau * v * v'
   /// where tau is a scalar, and v is a vector with
   ///
   ///         v[0] through v[i-1] = 0; v[i]   = 1
   /// with v[i+1] through v[m-1] stored on exit below the diagonal
   /// in the ith column of A, and tau in tau[i].
   /// @tparam real_t Floating point type.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @param m The number of rows of the matrix A.
   /// @param n The number of columns of the matrix A.
   /// @param[in,out] A Real m-by-n matrix.
   /// On exit, the elements on and above the diagonal of the array
   /// contain the min(m,n)-by-n upper trapezoidal matrix R
   /// (R is upper triangular if m >= n); the elements below the diagonal,
   /// with the array tau, represent the unitary matrix Q as a
   /// product of elementary reflectors.
   /// @param ldA Column length of the matrix A.  ldA>=m.
   /// @param[out] tau Real vector of length min(m,n).
   /// The scalar factors of the elementary reflectors.
   /// @param w Workspace vector of length n (optional).
   /// If not used, workspace is managed internally.

   template<typename real_t>
   int GRQR2(int_t m, int_t n, real_t *A, int_t ldA, real_t *tau, real_t *w=NULL)
   {
      using std::min;
      const real_t one(1.0);
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;

      int_t k=min(m,n);

      bool allocate=(w==NULL)?1:0;
      if(allocate)
         w=new real_t[n];
      real_t *B=A;
      for(int_t i=0;i<k;i++)
      {
         B+=ldA;
         LARFG(m-i,A[i],&(A[i+1]),1,tau[i]);
         if(i<n-1)
         {
            real_t alpha=A[i];
            A[i]=one;
            LARF('L',m-i,n-i-1,A[i],1,tau[i],B+i,ldA,w);
            A[i]=alpha;
         }
      }
      if(allocate)
         delete [] w;
      return 0;
   }

   /// @brief Computes a QR factorization of a complex matrix A.
   ///
   /// The matrix Q is represented as a product of elementary reflectors
   ///
   ///         Q = H_1 H_2 . . . H_k, where k = min(m,n).
   /// Each H_i has the form
   ///
   ///         H_i = I - tau * v * v'
   /// where tau is a scalar, and v is a vector with
   ///
   ///         v[0] through v[i-1] = 0; v[i]   = 1
   /// with v[i+1] through v[m-1] stored on exit below the diagonal
   /// in the ith column of A, and tau in tau[i].
   /// @tparam real_t Floating point type.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @param m The number of rows of the matrix A.
   /// @param n The number of columns of the matrix A.
   /// @param[in,out] A Complex m-by-n matrix.
   /// On exit, the elements on and above the diagonal of the array
   /// contain the min(m,n)-by-n upper trapezoidal matrix R
   /// (R is upper triangular if m >= n); the elements below the diagonal,
   /// with the array tau, represent the unitary matrix Q as a
   /// product of elementary reflectors.
   /// @param ldA Column length of the matrix A.  ldA>=m.
   /// @param[out] tau Complex vector of length min(m,n).
   /// The scalar factors of the elementary reflectors.
   /// @param w Workspace vector of length n (optional).
   /// If not used, workspace is managed internally.

   template<typename real_t>
   int GEQR2(int_t m, int_t n, complex<real_t> *A, int_t ldA, complex<real_t> *tau, complex<real_t> *w=NULL)
   {
      return GEQR2< complex<real_t> >(m,n,A,ldA,tau,w);
   }
}
#endif
