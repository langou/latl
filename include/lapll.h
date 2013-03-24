//
//  lapll.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/26/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lapll_h
#define _lapll_h

/// @file lapll.h Computes smallest singular value of a matrix formed by two vectors.

#include "axpy.h"
#include "larfg.h"
#include "las2.h"
#include "dot.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes smallest singular value of a matrix formed by two real vectors.
   ///
   /// Given two column vectors X and Y, let
   ///
   ///        A = ( x y )
   /// The routine first computes the QR factorization of A = Q*R,
   /// and then computes the SVD of the 2-by-2 upper triangular matrix R.
   /// The smaller singular value of R is returned in ssmin, which is used
   /// as the measurement of the linear dependency of the vectors x and y.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors x and y.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to real vector y.
   /// @param incy Increment of vector y.
   /// @param ssmin On exit, the smallest singular value of the matrix A.
   /// @ingroup VEC
   
   template<typename real_t>
   void LAPLL(int_t n,real_t *x,int_t incx,real_t *y,int_t incy,real_t &ssmin)
   {
      const real_t one=1.0;
      const real_t zero=0.0;
      ssmin=zero;
      if(n>0)
      {
         real_t tau;
         LARFG(n,x[0],x+incx,incx,tau);
         real_t a11=x[0];
         x[0]=one;
         real_t c=-tau*DOT(n,x,incx,y,incy);
         AXPY(n,c,x,incx,y,incy);
         LARFG(n-1,y[incy],y+2*incy,incy,tau);
         real_t a12=y[0];
         real_t a22=y[incy];
         real_t ssmax;
         LAS2(a11,a12,a22,ssmin,ssmax);
      }
   }

   /// @brief Computes smallest singular value of a matrix formed by two complex vectors.
   ///
   /// Given two column vectors X and Y, let
   ///
   ///        A = ( x y )
   /// The routine first computes the QR factorization of A = Q*R,
   /// and then computes the SVD of the 2-by-2 upper triangular matrix R.
   /// The smaller singular value of R is returned in ssmin, which is used
   /// as the measurement of the linear dependency of the vectors x and y.
   /// @tparam real_t Floating point type.
   /// @param n Length of vectors x and y.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of vector x.
   /// @param y Pointer to complex vector y.
   /// @param incy Increment of vector y.
   /// @param ssmin On exit, the smallest singular value of the matrix A.
   /// @ingroup VEC
   
   template<typename real_t>
   void LAPLL(int_t n,complex<real_t> *x,int_t incx,complex<real_t> *y,int_t incy,real_t &ssmin)
   {
      using std::conj;
      const real_t one=1.0;
      const real_t zero=0.0;
      ssmin=zero;
      if(n>0)
      {
         complex<real_t> tau;
         LARFG(n,x[0],x+incx,incx,tau);
         complex<real_t> a11=x[0];
         x[0]=one;
         complex<real_t> c=-conj(tau)*DOTC(n,x,incx,y,incy);
         AXPY(n,c,x,incx,y,incy);
         LARFG(n-1,y[incy],y+2*incy,incy,tau);
         complex<real_t> a12=y[0];
         complex<real_t> a22=y[incy];
         real_t ssmax;
         LAS2(abs(a11),abs(a12),abs(a22),ssmin,ssmax);
      }
   }
}

#endif
