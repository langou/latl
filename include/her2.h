//
//  her2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/23/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _her2_h
#define _her2_h

/// @file her2.h Performs complex outer product of two vectors.


#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs a vector outer product of two complex vectors.
   ///
   /// For a complex Hermitian matrix A, complex vectors x and y, and complex scalar alpha,
   ///
   ///        A := alpha * x * y' + conj(alpha) * y * x' + A
   ///
   /// is computed. The result is Hermitian and is stored as either upper or lower triangular in A.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   ///
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param alpha Complex scalar.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of the vector x.  x!=0
   /// @param y Pointer to complex vector y.
   /// @param incy Increment of the vector y.  y!=0
   /// @param A Pointer to complex Hermitian n-by-n matrix A.If uplo = 'U' or 'u' then only the upper
   /// triangular part of A is referenced and the lower part is not referenced.  If uplo = 'L' or 'l' then
   /// only the lower triangular part of A is referenced and the upper part is not referenced.
   /// @param ldA Column length of matrix A.  ldA>=n.
   /// @ingroup BLAS

   template <typename real_t>
   int HER2(char uplo, int_t n, complex<real_t> alpha, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy, complex<real_t> *A, int_t ldA)
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
      else if(ldA<n)
         return -9;
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
               A+=ldA;
            }
         }
         else
         {
            for(j=0;j<n;j++)
            {
               tx=conj(alpha*x[j]);
               ty=alpha*conj(y[j]);
               A[j]=real(A[j])+real(x[j]*ty+y[j]*tx);
               for(i=j+1;i<n;i++)
                  A[i]+=x[i]*ty+y[i]*tx;
               A+=ldA;
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
               A+=ldA;
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
               A[j]=real(A[j])+real(x[jx]*ty+y[jy]*tx);
               for(i=j+1;i<n;i++)
               {
                  ix+=incx;
                  iy+=incy;
                  A[i]+=x[ix]*ty+y[iy]*tx;
               }
               jx+=incx;
               jy+=incy;
               A+=ldA;
            }
         }
      }
      return 0;
   }
#ifdef __latl_cblas
#include <cblas.h>

   template <> int HER2<float>(char uplo, int_t n, complex<float> alpha, complex<float> *x, int_t incx, complex<float> *y, int_t incy, complex<float> *A, int_t ldA)
   {
      uplo=std::toupper(uplo);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(incx==0)
         return -5;
      else if(ldA<n)
         return -7;

      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;

      cblas_cher2(CblasColMajor,Uplo,n,&alpha,x,incx,y,incy,A,ldA);

      return 0;
   }

   template <> int HER2<double>(char uplo, int_t n, complex<double> alpha, complex<double> *x, int_t incx, complex<double> *y, int_t incy, complex<double> *A, int_t ldA)
   {
      uplo=std::toupper(uplo);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(incx==0)
         return -5;
      else if(ldA<n)
         return -7;

      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;

      cblas_zher2(CblasColMajor,Uplo,n,&alpha,x,incx,y,incy,A,ldA);
      
      return 0;
   }
#endif

}
#endif
