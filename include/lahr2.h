//
//  lahr2.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 2/28/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _lahr2_h
#define _lahr2_h

/// @file lahr2.h Reduces the specified number of first columns of a general 
/// rectangular matrix A so that elements below the specified subdiagonal are 
/// zero, and returns auxiliary matrices which are needed to apply the 
/// transformation to the unreduced part of A.

#include <algorithm>
#include "latl.h"
#include "gemv.h"
#include "copy.h"
#include "trmv.h"
#include "axpy.h"
#include "larfg.h"
#include "scal.h"
#include "lacpy.h"
#include "trmm.h"
#include "gemm.h"
#include "lacgv.h"

namespace LATL
{
   /// @brief Reduces the first NB columns of A real general n-BY-(n-k+1)
   /// matrix A so that elements below the k-th subdiagonal are zero. The 
   /// reduction is performed by an orthogonal similarity transformation 
   /// Q**T * A * Q. The routine returns the matrices V and T which determine Q
   /// as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.
   ///
   /// This is an auxiliary routine called by GEHRD.
   ///
   /// @return 0 if success.
   /// @param n Order of the matrix A.
   /// @param k The offset for the reduction. Elements below the k-th 
   /// subdiagonal in the first NB columns are reduced to zero.  k < n.
   /// @param nb The number of columns to be reduced.
   /// @param A Pointer to a real array, dimension (LDA,N-K+1). On entry, the
   /// n-by-(n-k+1) general matrix A.  On exit, the elements on and above the 
   /// k-th subdiagonal in the first NB columns are overwritten with the 
   /// corresponding elements of the reduced matrix; the elements below the 
   /// k-th subdiagonal, with the array TAU, represent the matrix Q as a 
   /// product of elementary reflectors. The other columns of A are unchanged.
   /// @param ldA Column length of matrix A.  ldA>=n
   /// @param tau Pointer to an array of dimension (nb). The scalar factors of
   /// the elementary reflector.
   /// @param T Pointer to a real array, dimension (ldT,nb) The upper
   /// triangular matrix T.
   /// @param ldT Column length of matrix T.  ldT>=nb
   /// @param Y Pointer to a real array, dimension (ldY,nb) 
   /// @param ldY Column length of matrix Y.  ldY>=n
   ///
   /// The matrix Q is represented as a product of nb elementary reflectors
   ///
   ///    Q = H(1) H(2) . . . H(nb).
   ///
   /// Each H(i) has the form
   ///
   ///    H(i) = I - tau * v * v**T
   ///
   /// where tau is a real scalar, and v is a real vector with
   /// v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
   /// A(i+k+1:n,i), and tau in TAU(i).
   ///
   /// The elements of the vectors v together form the (n-k+1)-by-nb matrix
   /// V which is needed, with T and Y, to apply the transformation to the
   /// unreduced part of the matrix, using an update of the form:
   /// A := (I - V*T*V**T) * (A - Y*V**T).
   ///
   /// The contents of A on exit are illustrated by the following example
   /// with n = 7, k = 3 and nb = 2:
   ///
   ///    ( a   a   a   a   a )
   ///    ( a   a   a   a   a )
   ///    ( a   a   a   a   a )
   ///    ( h   h   a   a   a )
   ///    ( v1  h   a   a   a )
   ///    ( v1  v2  a   a   a )
   ///    ( v1  v2  a   a   a )
   ///
   /// where a denotes an element of the original matrix A, h denotes a
   /// modified element of the upper Hessenberg matrix H, and vi denotes an
   /// element of the vector defining H(i).
   ///
   /// This subroutine is a slight modification of LAPACK-3.0's LAHRD
   /// incorporating improvements proposed by Quintana-Orti and Van de
   /// Gejin. Note that the entries of A(1:K,2:NB) differ from those
   /// returned by the original LAPACK-3.0's LAHRD routine. (This
   /// subroutine is not backward compatible with LAPACK-3.0's LAHRD.)
   ///
   ///References:
   ///================
   ///
   /// Gregorio Quintana-Orti and Robert van de Geijn, "Improving the
   /// performance of reduction to Hessenberg form," ACM Transactions on
   /// Mathematical Software, 32(2):180-194, June 2006.
   ///
   template <typename real_t>
   int_t lahr2( int_t n, int_t k, int_t nb, real_t *A, int_t ldA, real_t *tau, real_t *T, int_t ldT, real_t *Y, int_t ldY )
   {
      using std::min;
      int_t i;
      real_t ei;

      for(i=0;i<nb;i++)
      {
         if(i>0)
         {
            gemv<real_t>( 'N', n-k, i, -1.0, Y+k, ldY, A+k+i-1, ldA, 1.0, A+k+i*ldA, 1);
            copy<real_t>( i, A+k+i*ldA, 1, T+(nb-1)*ldT, 1);
            trmv<real_t>( 'L', 'T', 'U', i, A+k, ldA, T+(nb-1)*ldT, 1);
            gemv<real_t>( 'T', n-k-i, i, 1.0, A+k+i, ldA, A+k+i+i*ldA, 1, 1.0, T+(nb-1)*ldT, 1);
            trmv<real_t>( 'U', 'T', 'N', i, T, ldT, T+(nb-1)*ldT, 1);
            gemv<real_t>( 'N', n-k-i, i, -1.0, A+k+i, ldA, T+(nb-1)*ldT, 1, 1.0, A+k+i+i*ldA, 1);
            trmv<real_t>( 'L', 'N', 'U', i, A+k, ldA, T+(nb-1)*ldT, 1);
            axpy<real_t>( i, -1.0, T+(nb-1)*ldT, 1, A+k+i*ldA, 1);
            A[k+i-1+(i-1)*ldA] = ei;
         }
         larfg<real_t>( n-k-i, A[k+i+i*ldA], A+min(k+i+1,n-1)+i*ldA, 1, tau[i]);
         ei = A[k+i+i*ldA];
         A[k+i+i*ldA] = 1.0;
         gemv<real_t>( 'N', n-k, n-k-i, 1.0, A+k+(i+1)*ldA, ldA, A+k+i+i*ldA, 1, 0.0, Y+k+i*ldY, 1); 
         gemv<real_t>( 'T', n-k-i, i, 1.0, A+k+i, ldA, A+k+i+i*ldA, 1, 0.0, T+i*ldT, 1);
         gemv<real_t>( 'N', n-k, i, -1.0, Y+k, ldY, T+i*ldT, 1, 1.0, Y+k+i*ldY, 1);
         scal<real_t>( n-k, tau[i], Y+k+i*ldY, 1);
         scal<real_t>( i, -tau[i], T+i*ldT, 1);
         trmv<real_t>( 'U', 'N', 'N', i, T, ldT, T+i*ldT, 1);
         T[i+i*ldT] = tau[i];
      }
      A[k+nb-1+(nb-1)*ldA] = ei;
      lacpy<real_t>( 'A', k, nb, A+ldA, ldA, Y, ldY);
      trmm<real_t>( 'R', 'L', 'N', 'U', k, nb, 1.0, A+k, ldA, Y, ldY);
      if(n>k+nb)
         gemm<real_t>( 'N', 'N', k, nb, n-k-nb, 1.0, A+(nb+1)*ldA, ldA, A+k+nb, ldA, 1.0, Y, ldY);
      trmm<real_t>( 'R', 'U', 'N', 'N', k, nb, 1.0, T, ldT, Y, ldY);
      return 0;
   }

   /// @brief Reduces the first NB columns of A real general n-BY-(n-k+1)
   /// matrix A so that elements below the k-th subdiagonal are zero. The 
   /// reduction is performed by an orthogonal similarity transformation 
   /// Q**T * A * Q. The routine returns the matrices V and T which determine Q
   /// as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.
   ///
   /// This is an auxiliary routine called by GEHRD.
   ///
   /// @return 0 if success.
   /// @param n Order of the matrix A.
   /// @param k The offset for the reduction. Elements below the k-th 
   /// subdiagonal in the first NB columns are reduced to zero.  k < n.
   /// @param nb The number of columns to be reduced.
   /// @param A Pointer to a real array, dimension (LDA,N-K+1). On entry, the
   /// n-by-(n-k+1) general matrix A.  On exit, the elements on and above the 
   /// k-th subdiagonal in the first NB columns are overwritten with the 
   /// corresponding elements of the reduced matrix; the elements below the 
   /// k-th subdiagonal, with the array TAU, represent the matrix Q as a 
   /// product of elementary reflectors. The other columns of A are unchanged.
   /// @param ldA Column length of matrix A.  ldA>=n
   /// @param tau Pointer to an array of dimension (nb). The scalar factors of
   /// the elementary reflector.
   /// @param T Pointer to a real array, dimension (ldT,nb) The upper
   /// triangular matrix T.
   /// @param ldT Column length of matrix T.  ldT>=nb
   /// @param Y Pointer to a real array, dimension (ldY,nb) 
   /// @param ldY Column length of matrix Y.  ldY>=n
   ///
   /// The matrix Q is represented as a product of nb elementary reflectors
   ///
   ///    Q = H(1) H(2) . . . H(nb).
   ///
   /// Each H(i) has the form
   ///
   ///    H(i) = I - tau * v * v**T
   ///
   /// where tau is a real scalar, and v is a real vector with
   /// v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
   /// A(i+k+1:n,i), and tau in TAU(i).
   ///
   /// The elements of the vectors v together form the (n-k+1)-by-nb matrix
   /// V which is needed, with T and Y, to apply the transformation to the
   /// unreduced part of the matrix, using an update of the form:
   /// A := (I - V*T*V**T) * (A - Y*V**T).
   ///
   /// The contents of A on exit are illustrated by the following example
   /// with n = 7, k = 3 and nb = 2:
   ///
   ///    ( a   a   a   a   a )
   ///    ( a   a   a   a   a )
   ///    ( a   a   a   a   a )
   ///    ( h   h   a   a   a )
   ///    ( v1  h   a   a   a )
   ///    ( v1  v2  a   a   a )
   ///    ( v1  v2  a   a   a )
   ///
   /// where a denotes an element of the original matrix A, h denotes a
   /// modified element of the upper Hessenberg matrix H, and vi denotes an
   /// element of the vector defining H(i).
   ///
   /// This subroutine is a slight modification of LAPACK-3.0's LAHRD
   /// incorporating improvements proposed by Quintana-Orti and Van de
   /// Gejin. Note that the entries of A(1:K,2:NB) differ from those
   /// returned by the original LAPACK-3.0's LAHRD routine. (This
   /// subroutine is not backward compatible with LAPACK-3.0's LAHRD.)
   ///
   ///References:
   ///================
   ///
   /// Gregorio Quintana-Orti and Robert van de Geijn, "Improving the
   /// performance of reduction to Hessenberg form," ACM Transactions on
   /// Mathematical Software, 32(2):180-194, June 2006.
   ///
   template <typename real_t>
   int_t lahr2( int_t n, int_t k, int_t nb, complex<real_t> *A, int_t ldA, complex<real_t> *tau, complex<real_t> *T, int_t ldT, complex<real_t> *Y, int_t ldY )
   {
      using std::min;
      int_t i;
      complex<real_t> ei;

      for(i=0;i<nb;i++)
      {
         if(i>0)
         {
            lacgv<real_t>( i, A+k+i-1, ldA);
            gemv<real_t>( 'N', n-k, i, -1.0, Y+k, ldY, A+k+i-1, ldA, 1.0, A+k+i*ldA, 1);
            lacgv<real_t>( i, A+k+i-1, ldA);
            copy<real_t>( i, A+k+i*ldA, 1, T+(nb-1)*ldT, 1);
            trmv<real_t>( 'L', 'C', 'U', i, A+k, ldA, T+(nb-1)*ldT, 1);
            gemv<real_t>( 'C', n-k-i, i, 1.0, A+k+i, ldA, A+k+i+i*ldA, 1, 1.0, T+(nb-1)*ldT, 1);
            trmv<real_t>( 'U', 'C', 'N', i, T, ldT, T+(nb-1)*ldT, 1);
            gemv<real_t>( 'N', n-k-i, i, -1.0, A+k+i, ldA, T+(nb-1)*ldT, 1, 1.0, A+k+i+i*ldA, 1);
            trmv<real_t>( 'L', 'N', 'U', i, A+k, ldA, T+(nb-1)*ldT, 1);
            axpy<real_t>( i, -1.0, T+(nb-1)*ldT, 1, A+k+i*ldA, 1);
            A[k+i-1+(i-1)*ldA] = ei;
         }
         larfg<real_t>( n-k-i, A[k+i+i*ldA], A+min(k+i+1,n-1)+i*ldA, 1, tau[i]);
         ei = A[k+i+i*ldA];
         A[k+i+i*ldA] = 1.0;
         gemv<real_t>( 'N', n-k, n-k-i, 1.0, A+k+(i+1)*ldA, ldA, A+k+i+i*ldA, 1, 0.0, Y+k+i*ldY, 1); 
         gemv<real_t>( 'C', n-k-i, i, 1.0, A+k+i, ldA, A+k+i+i*ldA, 1, 0.0, T+i*ldT, 1);
         gemv<real_t>( 'N', n-k, i, -1.0, Y+k, ldY, T+i*ldT, 1, 1.0, Y+k+i*ldY, 1);
         scal<real_t>( n-k, tau[i], Y+k+i*ldY, 1);
         scal<real_t>( i, -tau[i], T+i*ldT, 1);
         trmv<real_t>( 'U', 'N', 'N', i, T, ldT, T+i*ldT, 1);
         T[i+i*ldT] = tau[i];
      }
      A[k+nb-1+(nb-1)*ldA] = ei;
      lacpy<real_t>( 'A', k, nb, A+ldA, ldA, Y, ldY);
      trmm<real_t>( 'R', 'L', 'N', 'U', k, nb, 1.0, A+k, ldA, Y, ldY);
      if(n>k+nb)
         gemm<real_t>( 'N', 'N', k, nb, n-k-nb, 1.0, A+(nb+1)*ldA, ldA, A+k+nb, ldA, 1.0, Y, ldY);
      trmm<real_t>( 'R', 'U', 'N', 'N', k, nb, 1.0, T, ldT, Y, ldY);
      return 0;
   }
}
#endif
