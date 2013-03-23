//
//  gehrd.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 2/28/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _gehrd_h
#define _gehrd_h

/// @file gehrd.h Reduces a general matrix A to upper Hessenberg form H by
/// an orthogonal (unitary) similarity transformation:  Q**T * A * Q = H .

#include <algorithm>
#include "latl.h"
#include "lahr2.h"
#include "gemm.h"
#include "trmm.h"
#include "axpy.h"
#include "larfb.h"
#include "gehd2.h"

namespace LATL
{
   /// @brief Reduces a real general matrix A to upper Hessenberg form H by an
   /// orthogonal similarity transformation:  Q**T * A * Q = H .
   /// @tparam real_t Floating point type.
   /// @param n The order of the matrix A.  n >= 0.
   /// @param ilo Column index.
   /// @param ihi Column index. It is assumed that A is already upper 
   /// triangular in rows and columns 1:ILO-1 and IHI+1:N. ILO and IHI are 
   /// normally set by a previous call to DGEBAL; otherwise they should be set
   /// to 0 and n-1 respectively.  0 <= ilo <= ihi <= n-1, if n > 0; ilo=1 and
   /// ihi=0, if n=0.
   /// @param[in,out] A Pointer to real array, dimension (ldA,n) On entry, the n-by-n 
   /// general matrix to be reduced.  On exit, the upper triangle and the first
   /// subdiagonal of A are overwritten with the upper Hessenberg matrix H, and
   /// the elements below the first subdiagonal, with the array tau, represent 
   /// the orthogonal matrix Q as a product of elementary reflectors.
   /// @param ldA Leading dimension of the matrix A.  ldA >= n
   /// @param[out] tau Pointer to array of dimension (n-1).  The scalar factors of
   /// the elementary reflectors.  Elements 0:ilo-1 and ihi:n-2 of tau are set
   /// to zero.
   /// @param nb Block size.  0<nb<=n
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
   /// exit in A(i+2:ihi,i), and tau in tau(i).
   /// 
   /// The contents of A are illustrated by the following example, with
   /// n = 7, ilo = 1 and ihi = 5:
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
   /// 
   /// This file is a slight modification of LAPACK-3.0's DGEHRD
   /// subroutine incorporating improvements proposed by Quintana-Orti and
   /// Van de Geijn (2006). (See DLAHR2.)
   template <typename real_t>
   int_t gehrd( int_t n, int_t ilo, int_t ihi, real_t *A, int_t ldA, real_t *tau, int_t nb=32)
   {

      using std::min;
      int_t i, j, nh, ib;
      int_t nx = 128;

      if(n<0)
         return -1;
      else if((ilo<0) || (ilo>n-1))
         return -2;
      else if((ihi<ilo) || (ihi>n-1))
         return -3;
      else if(ldA<n)
         return -5;

      real_t ei;
      real_t *T=new real_t[n*(nb+1)];
      real_t *work=new real_t[n*nb];

      for(i=0;i<ilo-1;i++)
         tau[i] = 0.0;
      for(i=ihi;i<n-1;i++)
         tau[i] = 0.0;

      nh = ihi - ilo + 1;

      if (nb >= nh)
         i = ilo;
      else
      {
         for(i=ilo;i<ihi-nx;i+=nb)
         {
            ib = min( nb, ihi-i );
            lahr2( ihi+1, i+1, ib, A+i*ldA, ldA, tau+i, T, n, work, n);
            ei = A[i+ib+(i+ib-1)*ldA];
            A[i+ib+(i+ib-1)*ldA] = 1.0;
            gemm<real_t>( 'N', 'T', ihi+1, ihi-i-ib+1, ib, -1.0, work, n, A+i+ib+i*ldA, ldA, 1.0, A+(i+ib)*ldA, ldA);
            A[i+ib+(i+ib-1)*ldA] = ei;
            trmm<real_t>( 'R', 'L', 'T', 'U', i+1, ib-1, 1.0, A+i+1+i*ldA, ldA, work, n);
            for(j=0;j<ib-1;j++)
               axpy<real_t>( i+1, -1.0, work+n*j, 1, A+(i+j+1)*ldA, 1);
            larfb( 'L', 'T', 'F', 'C', ihi-i, n-i-ib, ib, A+i+1+i*ldA, ldA, T, n, A+i+1+(i+ib)*ldA, ldA, work);
         }
      }

      gehd2( n, i, ihi, A, ldA, tau );

      delete [] work;
      return 0;

   }

   /// @brief Reduces a complex general matrix A to upper Hessenberg form H by an
   /// unitary similarity transformation:  Q**T * A * Q = H .
   /// @tparam real_t Floating point type.
   /// @param n The order of the matrix A.  n >= 0.
   /// @param ilo Column index.
   /// @param ihi Column index. It is assumed that A is already upper 
   /// triangular in rows and columns 0:ilo-1 and ihi+1:n-1. ilo and ihi are 
   /// normally set by a previous call to DGEBAL; otherwise they should be set
   /// to 0 and n-1 respectively. 0 <= ilo <= ihi <= n-1, if n > 0; ilo=1 and 
   /// ihi=0, if n=0.
   /// @param[in,out] A Pointer to complex array, dimension (ldA,n) On entry, 
   /// the n-by-n general matrix to be reduced.  On exit, the upper triangle
   /// and the first subdiagonal of A are overwritten with the upper Hessenberg
   /// matrix H, and the elements below the first subdiagonal, with the array
   /// tau, represent the unitary matrix Q as a product of elementary 
   /// reflectors.
   /// @param ldA Leading dimension of the matrix A.  ldA >= n
   /// @param[out] tau Pointer to array of dimension (n-1).  The scalar factors of
   /// the elementary reflectors.  Elements 0:ilo-1 and ihi:n-2 of tau are set
   /// to zero.
   /// @param nb Block size.  0<nb<=n
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
   /// exit in A(i+2:ihi,i), and tau in tau(i).
   /// 
   /// The contents of A are illustrated by the following example, with
   /// n = 7, ilo = 1 and ihi = 5:
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
   /// 
   /// This file is a slight modification of LAPACK-3.0's DGEHRD
   /// subroutine incorporating improvements proposed by Quintana-Orti and
   /// Van de Geijn (2006). (See DLAHR2.)
   template <typename real_t>
   int_t gehrd( int_t n, int_t ilo, int_t ihi, complex<real_t> *A, int_t ldA, complex<real_t> *tau, int_t nb=32)
   {

      using std::min;
      int_t i, j, nh, ib;
      int_t nx = 128;

      if(n<0)
         return -1;
      else if((ilo<0) || (ilo>n-1))
         return -2;
      else if((ihi<ilo) || (ihi>n-1))
         return -3;
      else if(ldA<n)
         return -5;

      complex<real_t> ei;
      complex<real_t> *T=new complex<real_t>[n*(nb+1)];
      complex<real_t> *work=new complex<real_t>[n*nb];

      for(i=0;i<ilo-1;i++)
         tau[i] = 0.0;
      for(i=ihi;i<n-1;i++)
         tau[i] = 0.0;

      nh = ihi - ilo + 1;

      if (nb >= nh)
         i = ilo;
      else
      {
         for(i=ilo;i<ihi-nx;i+=nb)
         {
            ib = min( nb, ihi-i );
            lahr2( ihi+1, i+1, ib, A+i*ldA, ldA, tau+i, T, n, work, n);
            ei = A[i+ib+(i+ib-1)*ldA];
            A[i+ib+(i+ib-1)*ldA] = 1.0;
            gemm<real_t>( 'N', 'C', ihi+1, ihi-i-ib+1, ib, -1.0, work, n, A+i+ib+i*ldA, ldA, 1.0, A+(i+ib)*ldA, ldA);
            A[i+ib+(i+ib-1)*ldA] = ei;
            trmm<real_t>( 'R', 'L', 'C', 'U', i+1, ib-1, 1.0, A+i+1+i*ldA, ldA, work, n);
            for(j=0;j<ib-1;j++)
               axpy<real_t>( i+1, -1.0, work+n*j, 1, A+(i+j+1)*ldA, 1);
            larfb( 'L', 'C', 'F', 'C', ihi-i, n-i-ib, ib, A+i+1+i*ldA, ldA, T, n, A+i+1+(i+ib)*ldA, ldA, work);
         }
      }

      gehd2( n, i, ihi, A, ldA, tau );

      delete [] work;
      return 0;

   }
}
#endif
