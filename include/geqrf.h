//
//  geqrf.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/26/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _geqrf_h
#define _geqrf_h

/// @file geqrf.h Computes a blocked QR factorization of a matrix.

#include <algorithm>
#include "geqr2.h"
#include "larfb.h"
#include "larft.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes a blocked QR factorization of a real matrix A.
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
   /// @param nb Blocking factor.  nb>0.
   /// @ingroup COMP

   template<typename real_t>
   int GEQRF(int_t m, int_t n,real_t *A, int_t ldA, real_t *tau, int_t nb)
   {
      using std::min;
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;
      else if(nb<1)
         return -6;

      int_t k=min(m,n);
      nb=min(nb,k);
      
      if(k==0)
         return 0;
      
      real_t *w=new real_t[n*nb+nb];
      real_t *B=A;
      for(int_t i=0;i<k;i+=nb)
      {
         int_t ib=min(k-i,nb);
         B+=ib*ldA;
         GEQR2(m-i,ib,A+i,ldA,tau);
         if(i+ib<n-1)
         {
            LARFT('F','C',m-i,ib,A+i,ldA,tau,w,ib);
            LARFB('L','T','F','C',m-i,n-i-ib,ib,A+i,ldA,w,ib,B+i,ldA,w+ib*ib);
         }
         A=B;
         tau+=ib;
      }
      delete [] w;
      return 0;
   }

   /// @brief Computes a blocked QR factorization of a complex matrix A.
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
   /// @param nb Blocking factor.  nb>0.
   /// @ingroup COMP

   template<typename real_t>
   int GEQRF(int_t m, int_t n, complex<real_t> *A, int_t ldA, complex<real_t> *tau, int_t nb)
   {
      using std::min;
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;
      else if(nb<1)
         return -6;

      int_t k=min(m,n);
      nb=min(nb,k);

      if(k==0)
         return 0;

      complex<real_t> *w=new complex<real_t>[n*nb+nb];
      complex<real_t> *B=A;
      for(int_t i=0;i<k;i+=nb)
      {
         int_t ib=min(k-i,nb);
         B+=ib*ldA;
         GEQR2(m-i,ib,A+i,ldA,tau);
         if(i+ib<n-1)
         {
            LARFT('F','C',m-i,ib,A+i,ldA,tau,w,ib);
            LARFB('L','C','F','C',m-i,n-i-ib,ib,A+i,ldA,w,ib,B+i,ldA,w+ib*ib);
         }
         A=B;
         tau+=ib;
      }
      delete [] w;
      return 0;
   }
}
#endif

