//
//  hpmv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/13/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _hpmv_h
#define _hpmv_h

/// @file hpmv.h Performs complex Hermitian matrix-vector multiplication using packed storage.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs complex Hermitian matrix-vector multiplication using packed storage.
   /// 
   /// For complex Hermitian matrix A, complex vectors x and y, and complex scalars alpha and beta
   ///
   ///         y := alpha*A*x+beta*y
   /// is computed, where A uses packed storage format.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the packed matrix A is upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular,
   ///        if uplo = 'L' or 'l' then A is lower triangular.
   /// @param n Specifies the order of the Hermitian matrix A.  n>=0
   /// @param alpha Complex scalar.
   /// @param A Pointer to packed complex Hermitian n-by-n matrix A.  The packed storage format stores a triangular matrix
   /// using variable length columns so that no space is wasted.  If A is upper triangular, the first column consists only
   /// of the diagonal element, which has length 1.  The second column has length 2, the third length 3, and so on until the
   /// nth column, which has length n.  Similarly, if A is lower triangular, the first column has length n, the second length n-1,
   /// and so on until the nth column, which has length 1.  The entire n-by-n triangular matrix is stored using n(n+1)/2 elements.
   /// @param x Pointer the complex vector x.
   /// @param incx Increment of vector x. incx!=0
   /// @param beta Complex scalar.
   /// @param y Pointer the complex vector y.
   /// @param incy Increment of vector y. incy!=0
   /// @ingroup BLAS

   template<typename real_t>
   int HPMV(char uplo, int_t n, complex<real_t> alpha, complex<real_t> *A, complex<real_t> *x, int_t incx, complex<real_t> beta, complex<real_t> *y, int_t incy)
   {
      using std::conj;
      using std::real;
      using std::toupper;
      const complex<real_t> one(1.0,0.0);
      const complex<real_t> zero(0.0,0.0);
      int_t kx,ky,i,ix,iy,j,jx,jy;
      complex<real_t> temp,sum;
      
      uplo=toupper(uplo);

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(incx==0)
         return -6;
      else if(incy==0)
         return -9;
      else if(n==0)
         return 0;
      
      kx=(incx>0)?0:(1-n)*incx;
      ky=(incy>0)?0:(1-n)*incy;
      
      if(beta==zero)
      {
         if(incy==1)
         {
            for(i=0;i<n;i++)
               y[i]=zero;
         }
         else
         {
            iy=ky;
            for(i=0;i<n;i++)
            {
               y[iy]=zero;
               iy+=incy;
            }
            
         }
      }
      else if(beta!=one)
      {
         if(incy==1)
         {
            for(i=0;i<n;i++)
               y[i]*=beta;
         }
         else
         {
            iy=ky;
            for(i=0;i<n;i++)
            {
               y[iy]*=beta;
               iy+=incy;
            }
         }
      }
      
      if(alpha!=zero)
      {
         if(uplo=='U')
         {
            if((incx==1)&&(incy==1))
            {
               for(j=0;j<n;j++)
               {
                  temp=alpha*x[j];
                  sum=zero;
                  for(i=0;i<j;i++)
                  {
                     y[i]+=temp*A[i];
                     sum+=conj(A[i])*x[i];
                  }
                  y[j]+=temp*real(A[j])+alpha*sum;
                  A+=j+1;
               }
            }
            else
            {
               jx=kx;
               jy=ky;
               for(j=0;j<n;j++)
               {
                  temp=alpha*x[jx];
                  sum=zero;
                  ix=kx;
                  iy=ky;
                  for(i=0;i<j;i++)
                  {
                     y[iy]+=temp*A[i];
                     sum+=conj(A[i])*x[ix];
                     ix+=incx;
                     iy+=incy;
                  }
                  y[jy]+=temp*real(A[j])+alpha*sum;
                  A+=j+1;
                  jx+=incx;
                  jy+=incy;
               }
            }
         }
         else
         {
            if((incx==1)&&(incy==1))
            {
               for(j=0;j<n;j++)
               {
                  temp=alpha*x[j];
                  sum=zero;
                  y[j]+=temp*real(A[0]);
                  for(i=j+1;i<n;i++)
                  {
                     y[i]+=temp*A[i-j];
                     sum+=conj(A[i-j])*x[i];
                  }
                  y[j]+=alpha*sum;
                  A+=n-j;
               }
            }
            else
            {
               jx=kx;
               jy=ky;
               for(j=0;j<n;j++)
               {
                  temp=alpha*x[jx];
                  sum=zero;
                  ix=jx;
                  iy=jy;
                  y[jy]+=temp*real(A[0]);
                  for(i=j+1;i<n;i++)
                  {
                     ix+=incx;
                     iy+=incy;
                     y[iy]+=temp*A[i-j];
                     sum+=conj(A[i-j])*x[ix];
                  }
                  y[jy]+=alpha*sum;
                  jx+=incx;
                  jy+=incy;
                  A+=n-j;
               }
            }
         }      
      }
      return 0;
   }
}
#endif

