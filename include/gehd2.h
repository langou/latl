//
//  gehd2.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 2/27/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _gehd2_h
#define _gehd2_h

/// @file gehd2.h Reduces a general matrix A to upper Hessenberg form H by an orthogonal 
/// similarity transformation:  Q**T * A * Q = H

#include <algorithm>
#include "latl.h"
#include "larfg.h"
#include "larf.h"

namespace latl
{
   /// @brief Reduces a general matrix A to upper Hessenberg form H by an orthogonal 
   /// similarity transformation:  Q**T * A * Q = H
   ///
   /// @return 0 if success.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param ilo An integer.
   /// @param ihi An integer.
   /// @param A Pointer to general matrix.  On entry, the n by n general matrix 
   /// to be reduced.  On exit, the upper triangle and the first subdiagonal of 
   /// A are overwritten with the upper Hessenberg matrix H, and the elements 
   /// below the first subdiagonal, with the array TAU, represent the orthogonal 
   /// matrix Q as a product of elementary reflectors. 
   /// @param ldA Column length of matrix A.  ldA>=n
   /// @param tau Pointer to an array of dimension (N-1). The scalar factors of
   /// the elementary reflector.
   ///
   /// It is assumed that A is already upper triangular in rows
   /// and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
   /// set by a previous call to GEBAL; otherwise they should be
   /// set to 1 and N respectively. 
   /// 0 <= ILO <= IHI < N
   ///
   /// The matrix Q is represented as a product of (ihi-ilo) elementary
   /// reflectors
   ///
   ///    Q = H(ilo) H(ilo+1) . . . H(ihi-1).
   ///
   /// Each H(i) has the form
   ///
   ///    H(i) = I - tau * v * v**T
   ///
   /// where tau is a real scalar, and v is a real vector with
   /// v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
   /// exit in A(i+2:ihi,i), and tau in TAU(i).
   ///
   /// The contents of A are illustrated by the following example, with
   /// n = 7, ilo = 2 and ihi = 6:
   ///
   /// on entry,                        on exit,
   ///
   /// ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
   /// (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
   /// (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
   /// (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
   /// (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
   /// (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
   /// (                         a )    (                          a )
   ///
   /// where a denotes an element of the original matrix A, h denotes a
   /// modified element of the upper Hessenberg matrix H, and vi denotes an
   /// element of the vector defining H(i).
   
   template <typename real_t>
   int_t gehd2(int_t n, int_t ilo, int_t ihi, real_t *A, int_t ldA, real_t *tau )
   {

      int_t i;
      real_t alpha;

      using std::min; 
      using latl::larfg;
      using latl::larf;

      if (n<0)
         return -1;
      else if((ilo<1) || (ilo>n))
         return -2;
      else if((ihi<ilo) || (ihi>n))
         return -3;
      else if(ldA<n)
         return -5;

      for(i=ilo;i<ihi;i++)
      {
         alpha = A[i+(i-1)*ldA];
         larfg<real_t>( ihi-i, alpha, A+min(i+2,n)-1+(i-1)*ldA, 1, tau[i-1]);
         A[i+(i-1)*ldA] = 1.0;

//          larf<real_t>( 'R', ihi, ihi-i, A+i+(i-1)*ldA, 1, tau[i-1], A+i*ldA, ldA);
         real_t *w1=new real_t[ihi];
         gemv<real_t>('N',ihi,ihi-i,1.0,A+i*ldA,ldA,A+i+(i-1)*ldA,1,0.0,w1,1);
         ger<real_t>(ihi,ihi-i,-tau[i-1],w1,1,A+i+(i-1)*ldA,1,A+i*ldA,ldA);
         delete [] w1;

//          larf<real_t>( 'L', ihi-i, n-i, A+i+(i-1)*ldA, 1, tau[i-1], A+i+i*ldA, ldA);
         real_t *w2=new real_t[n-i];
         gemv<real_t>('T',ihi-i,n-i,1.0,A+i+i*ldA,ldA,A+i+(i-1)*ldA,1,0.0,w2,1);
         ger<real_t>(ihi-i,n-i,-tau[i-1],A+i+(i-1)*ldA,1,w2,1,A+i+i*ldA,ldA);
         delete [] w2;
         A[i+(i-1)*ldA] = alpha;
      }
      return 0;
   }

   /// @brief Reduces a general matrix A to upper Hessenberg form H by an orthogonal 
   /// similarity transformation:  Q**T * A * Q = H
   ///
   /// @return 0 if success.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param ilo An integer.
   /// @param ihi An integer.
   /// @param A Pointer to general matrix.  On entry, the n by n general matrix 
   /// to be reduced.  On exit, the upper triangle and the first subdiagonal of 
   /// A are overwritten with the upper Hessenberg matrix H, and the elements 
   /// below the first subdiagonal, with the array TAU, represent the orthogonal 
   /// matrix Q as a product of elementary reflectors. 
   /// @param ldA Column length of matrix A.  ldA>=n
   /// @param tau Pointer to an array of dimension (N-1). The scalar factors of
   /// the elementary reflector.
   ///
   /// It is assumed that A is already upper triangular in rows
   /// and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
   /// set by a previous call to GEBAL; otherwise they should be
   /// set to 1 and N respectively. 
   /// 0 <= ILO <= IHI < N
   ///
   /// The matrix Q is represented as a product of (ihi-ilo) elementary
   /// reflectors
   ///
   ///    Q = H(ilo) H(ilo+1) . . . H(ihi-1).
   ///
   /// Each H(i) has the form
   ///
   ///    H(i) = I - tau * v * v**H
   ///
   /// where tau is a complex scalar, and v is a complex vector with
   /// v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
   /// exit in A(i+2:ihi,i), and tau in TAU(i).
   ///
   /// The contents of A are illustrated by the following example, with
   /// n = 7, ilo = 2 and ihi = 6:
   ///
   /// on entry,                        on exit,
   ///
   /// ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
   /// (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
   /// (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
   /// (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
   /// (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
   /// (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
   /// (                         a )    (                          a )
   ///
   /// where a denotes an element of the original matrix A, h denotes a
   /// modified element of the upper Hessenberg matrix H, and vi denotes an
   /// element of the vector defining H(i).
   
   template <typename real_t>
   int_t gehd2(int_t n, int_t ilo, int_t ihi, complex<real_t> *A, int_t ldA, complex<real_t> *tau )
   {

      int_t i;
      complex<real_t> alpha;

      using std::min; 
      using latl::larfg;
      using latl::larf;

      if (n<0)
         return -1;
      else if((ilo<1) || (ilo>n))
         return -2;
      else if((ihi<ilo) || (ihi>n))
         return -3;
      else if(ldA<n)
         return -5;

      for(i=ilo;i<ihi;i++)
      {
         alpha = A[i+(i-1)*ldA];
         larfg< real_t >( ihi-i, alpha, A+min(i+2,n)-1+(i-1)*ldA, 1, tau[i-1]);
         A[i+(i-1)*ldA] = 1.0;
         larf< real_t >( 'R', ihi, ihi-i, A+i+(i-1)*ldA, 1, tau[i-1], A+i*ldA, ldA);
         larf< real_t >( 'L', ihi-i, n-i, A+i+(i-1)*ldA, 1, std::conj(tau[i-1]), A+i+i*ldA, ldA);
         A[i+(i-1)*ldA] = alpha;
      }
      return 0;
   }
}
#endif
