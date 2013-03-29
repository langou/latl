//
//  hpr2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/23/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _hpr2_h
#define _hpr2_h

/// @file hpr2.h Performs complex outer products of two vectors using packed format.


#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs a complex vector outer product using packed storage.
   ///
   /// For complex vectors x and y, and complex scalar alpha
   ///
   ///         A := alpha*x*y'+conj(alpha)*y*x'+A,
   /// is computed, where the complex Hermitian matrix A uses packed storage format.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the packed matrix A is upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular,
   ///        if uplo = 'L' or 'l' then A is lower triangular.
   /// @param n Specifies the order of the Hermitian matrix A.  n>=0
   /// @param alpha Complex scalar.
   /// @param x Pointer the complex vector x.
   /// @param incx Increment of vector x. incx!=0
   /// @param y Pointer the complex vector y.
   /// @param incy Increment of vector y. incy!=0
   /// @param A Pointer to packed complex Hermitian n-by-n matrix A.  The packed storage format stores a triangular matrix
   /// using variable length columns so that no space is wasted.  If A is upper triangular, the first column consists only
   /// of the diagonal element, which has length 1.  The second column has length 2, the third length 3, and so on until the
   /// nth column, which has length n.  Similarly, if A is lower triangular, the first column has length n, the second length n-1,
   /// and so on until the nth column, which has length 1.  The entire n-by-n triangular matrix is stored using n(n+1)/2 elements.
   /// @ingroup BLAS

   template <typename real_t>
   int HPR2(char uplo, int_t n, complex<real_t> alpha, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy, complex<real_t> *A)
   {
      using std::conj;
      using std::real;
      using std::toupper;
      const complex<real_t> zero(0.0,0.0);
      complex<real_t> tx,ty;
      int_t i,j,kx,ky,jx,ix,jy,iy;

      uplo=toupper(uplo);

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(incx==0)
         return -5;
      else if(incy==0)
         return -7;
      else if((n==0)||(alpha==zero))
         return 0;

      if((incx==1)&&(incy==1))
      {
         if(uplo=='U')
         {
            for(j=0;j<n;j++)
            {
               tx=conj(alpha*x[j]);
               ty=alpha*conj(y[j]);
               for(i=0;i<j;i++)
                  A[i]+=x[i]*ty+y[i]*tx;
               A[j]=real(A[j])+real(x[j]*ty+y[j]*tx);
               A+=j+1;
            }
         }
         else
         {
            for(j=0;j<n;j++)
            {
               tx=conj(alpha*x[j]);
               ty=alpha*conj(y[j]);
               A[0]=real(A[0])+real(x[j]*ty+y[j]*tx);
               for(i=j+1;i<n;i++)
                  A[i-j]+=x[i]*ty+y[i]*tx;
               A+=n-j;
            }
         }
      }
      else
      {
         kx=(incx>0)?0:(1-n)*incx;
         ky=(incy>0)?0:(1-n)*incy;
         if(uplo=='U')
         {
            jx=kx;
            jy=ky;
            for(j=0;j<n;j++)
            {
               tx=conj(alpha*x[jx]);
               ty=alpha*conj(y[jy]);
               ix=kx;
               iy=ky;
               for(i=0;i<j;i++)
               {
                  A[i]+=x[ix]*ty+y[iy]*tx;
                  ix+=incx;
                  iy+=incy;
               }
               A[j]=real(A[j])+real(x[jx]*ty+y[jy]*tx);
               A+=j+1;
               jx+=incx;
               jy+=incy;
            }
         }
         else
         {
            jx=kx;
            jy=ky;
            for(j=0;j<n;j++)
            {
               tx=conj(alpha*x[jx]);
               ty=alpha*conj(y[jy]);
               ix=jx;
               iy=jy;
               A[0]=real(A[0])+real(x[jx]*ty+y[jy]*tx);
               for(i=j+1;i<n;i++)
               {
                  ix+=incx;
                  iy+=incy;
                  A[i-j]+=x[ix]*ty+y[iy]*tx;
               }
               jx+=incx;
               jy+=incy;
               A+=n-j;
            }
         }
      }
      return 0;
   }
}
#endif
