//
//  syr.h
//  Linear Algebra Template Library  
//
//  Created by Rodney James on 12/23/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _syr_h
#define _syr_h

/// @file syr.h Performs vector outer product.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs a vector outer product of a real vector with itself.
   /// 
   /// For a real symmetric matrix A, real vector x and real scalar alpha,
   ///
   ///        A := alpha * x * x' + A
   ///
   /// is computed.  The result is symmetric and is stored as either upper or lower triangular in A.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   ///
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param alpha Real scalar.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of the vector x.  x!=0
   /// @param A Pointer to real symmetric n-by-n matrix A.  If uplo = 'U' or 'u' then only the upper
   /// triangular part of A is referenced and the lower part is not referenced.  If uplo = 'L' or 'l' then
   /// only the lower triangular part of A is referenced and the upper part is not referenced.
   /// @param ldA Column length of matrix A.  ldA>=n.
   /// @ingroup VEC

   template <typename real_t>
   int syr(char uplo, int_t n, real_t alpha, real_t *x, int_t incx, real_t *A, int_t ldA)
   {
      using std::toupper;

      const real_t zero(0.0);
      real_t t;
      int_t i,j,kx,jx,ix;

      uplo=toupper(uplo);
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(incx==0)
         return -5;
      else if(ldA<n)
         return -7;
      else if((n==0)||(alpha==zero))
         return 0;
      
      if(incx==1)
      {
         if(uplo=='U')
         {
            for(j=0;j<n;j++)
            {
               t=alpha*x[j];
               for(i=0;i<=j;i++)
                  A[i]+=x[i]*t;
               A+=ldA;
            }
         }
         else
         {
            for(j=0;j<n;j++)
            {
               t=alpha*x[j];
               for(i=j;i<n;i++)
                  A[i]+=x[i]*t;
               A+=ldA;
            }
         }
      }
      else
      {
         kx=(incx>0)?0:(1-n)*incx;
         if(uplo=='U')
         {
            jx=kx;
            for(j=0;j<n;j++)
            {
               t=alpha*x[jx];
               ix=kx;
               for(i=0;i<=j;i++)
               {
                  A[i]+=x[ix]*t;
                  ix+=incx;
               }
               A+=ldA;
               jx+=incx;
            }
         }
         else
         {
            jx=kx;
            for(j=0;j<n;j++)
            {
               t=alpha*x[jx];
               ix=jx;
               for(i=j;i<n;i++)
               {
                  A[i]+=x[ix]*t;
                  ix+=incx;
               }
               jx+=incx;
               A+=ldA;
            }
         }
      }
      return 0;
   }
   
   /// @brief Performs a vector outer product of a complex vector with itself.
   /// 
   /// For a complex symmetric matrix A, complex vector x and complex scalar alpha,
   ///
   ///        A := alpha * x * x.' + A
   ///
   /// is computed.  The result is symmetric and is stored as either upper or lower triangular in A.
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
   /// @param A Pointer to complex symmetric n-by-n matrix A.  If uplo = 'U' or 'u' then only the upper
   /// triangular part of A is referenced and the lower part is not referenced.  If uplo = 'L' or 'l' then
   /// only the lower triangular part of A is referenced and the upper part is not referenced.
   /// @param ldA Column length of matrix A.  ldA>=n.
   /// @ingroup VEC

   template <typename real_t>
   int syr(char uplo, int_t n, complex<real_t> alpha, complex<real_t> *x, int_t incx, complex<real_t> *A, int_t ldA)
   {
      return syr< complex<real_t> >(uplo,n,alpha,x,incx,A,ldA);
   }
   
   /// @brief Performs a vector outer product of two real vectors.  
   /// 
   /// For a real symmetric matrix A, real vectors x and y, and real scalar alpha,
   ///
   ///        A := alpha * x * y' + alpha * y * x' + A
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
   /// @param alpha Real scalar.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of the vector x.  x!=0
   /// @param y Pointer to real vector y.
   /// @param incy Increment of the vector y.  y!=0
   /// @param A Pointer to real symmetric n-by-n matrix A.If uplo = 'U' or 'u' then only the upper
   /// triangular part of A is referenced and the lower part is not referenced.  If uplo = 'L' or 'l' then
   /// only the lower triangular part of A is referenced and the upper part is not referenced.
   /// @param ldA Column length of matrix A.  ldA>=n.
   /// @ingroup VEC

   template <typename real_t>
   int syr(char uplo, int_t n, real_t alpha, real_t *x, int_t incx, real_t *y, int_t incy, real_t *A, int_t ldA)
   {
      using std::toupper;

      const real_t zero(0.0);
      real_t tx,ty;
      int_t  i,j,kx,ky,jx,ix,jy,iy;

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
               tx=alpha*x[j];
               ty=alpha*y[j];
               for(i=0;i<=j;i++)
                  A[i]+=x[i]*ty+y[i]*tx;
               A+=ldA;
            }
         }
         else
         {
            for(j=0;j<n;j++)
            {
               tx=alpha*x[j];
               ty=alpha*y[j];
               for(i=j;i<n;i++)
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
               tx=alpha*x[jx];
               ty=alpha*y[jy];
               ix=kx;
               iy=ky;
               for(i=0;i<=j;i++)
               {
                  A[i]+=x[ix]*ty+y[iy]*tx;
                  ix+=incx;
                  iy+=incy;
               }
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
               tx=alpha*x[jx];
               ty=alpha*y[jy];
               ix=jx;
               iy=jy;
               for(i=j;i<n;i++)
               {
                  A[i]+=x[ix]*ty+y[iy]*tx;
                  ix+=incx;
                  iy+=incy;
               }
               jx+=incx;
               jy+=incy;
               A+=ldA;
            }
         }
      }
      return 0;
   }
}
#endif
