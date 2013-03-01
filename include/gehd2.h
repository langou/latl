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
#include "gemv.h"
#include "ger.h"

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
   /// set to 0 and N-1 respectively. 
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

      using std::min; 
      using std::max; 
      using latl::larfg;

      if (n<0)
         return -1;
      else if((ilo<0) || (ilo>n-1))
         return -2;
      else if((ihi<ilo) || (ihi>n-1))
         return -3;
      else if(ldA<n)
         return -5;

      int_t i;
      real_t alpha;
      real_t *w=new real_t[max(n-ilo-1,ihi+1)];
      real_t *v = A+ilo+1+ilo*ldA;
      real_t *CL = A+(ilo+1)*ldA;
      real_t *CR = v+ldA;
      for(i=ilo;i<ihi;i++)
      {
         alpha = v[0];
         larfg<real_t>( ihi-i, alpha, A+min(i+2,n-1)+i*ldA, 1, tau[i]);
         v[0] = 1.0;
         gemv<real_t>('N',ihi+1,ihi-i,1.0,CL,ldA,v,1,0.0,w,1);
         ger<real_t>(ihi+1,ihi-i,-tau[i],w,1,v,1,CL,ldA);
         gemv<real_t>('T',ihi-i,n-i-1,1.0,CR,ldA,v,1,0.0,w,1);
         ger<real_t>(ihi-i,n-i-1,-tau[i],v,1,w,1,CR,ldA);
         v[0] = alpha;
         v += 1+ldA;
         CL += ldA;
         CR = v+ldA;
      }
      delete [] w;
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
   /// set to 0 and N-1 respectively. 
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

      using std::min; 
      using std::max; 
      using std::conj;
      using latl::larfg;

      if (n<0)
         return -1;
      else if((ilo<0) || (ilo>n-1))
         return -2;
      else if((ihi<ilo) || (ihi>n-1))
         return -3;
      else if(ldA<n)
         return -5;

      int_t i;
      complex<real_t> alpha;
      complex<real_t> *w=new complex<real_t>[max(n-ilo-1,ihi+1)];
      complex<real_t> *v = A+ilo+1+ilo*ldA;
      complex<real_t> *CL = A+(ilo+1)*ldA;
      complex<real_t> *CR = v+ldA;
      for(i=ilo;i<ihi;i++)
      {
         alpha = v[0];
         larfg< real_t >( ihi-i, alpha, A+min(i+2,n-1)+i*ldA, 1, tau[i]);
         v[0] = 1.0;
         gemv<real_t>('N',ihi+1,ihi-i,1.0,CL,ldA,v,1,0.0,w,1);
         gerc<real_t>(ihi+1,ihi-i,-tau[i],w,1,v,1,CL,ldA);
         gemv<real_t>('C',ihi-i,n-i-1,1.0,CR,ldA,v,1,0.0,w,1);
         gerc<real_t>(ihi-i,n-i-1,-conj(tau[i]),v,1,w,1,CR,ldA);
         v[0] = alpha;
         v += 1+ldA;
         CL += ldA;
         CR = v+ldA;
      }
      delete [] w;
      return 0;
   }
}
#endif
