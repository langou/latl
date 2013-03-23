//
//  symv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/12/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _symv_h
#define _symv_h

/// @file symv.h Performs symmetric matrix-vector multiplication.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs real symmetric matrix-vector multiplication.
   ///
   /// For real symmetric matrix A, real vectors x and y, and real scalars alpha and beta,
   ///
   ///         y := alpha*A*x + beta*y
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetric matrix A
   /// is to be referenced:  
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular,
   ///        if uplo = 'L' or 'l' then A is lower triangular.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param alpha Real scalar.
   /// @param A Pointer to real n-by-n symmetric matrix A.  
   /// Only the upper or lower triangular part of A is referenced, depending on the value of uplo above.
   /// @param ldA Column length of the matrix A.  ldA>=n
   /// @param x Pointer to real vector x.
   /// @param incx Increment of the vector x.
   /// @param beta Real scalar.
   /// @param y Pointer to real vector y.
   /// @param incy Increment of the vector y.
   /// @ingroup MATV

   template <typename real_t>
   int SYMV(char uplo, int_t n, real_t alpha, real_t *A, int_t ldA, real_t *x, int_t incx, real_t beta, real_t *y, int_t incy)
   {
      using std::toupper;

      const real_t one(1.0);
      const real_t zero(0.0);
      int_t kx,ky,i,ix,iy,j,jx,jy;
      real_t temp,sum;

      uplo=toupper(uplo);

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -5;
      else if(incx==0)
         return -7;
      else if(incy==0)
         return -10;
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
                     sum+=A[i]*x[i];
                  }
                  y[j]+=temp*A[j]+alpha*sum;
                  A+=ldA;
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
                     sum+=A[i]*x[ix];
                     ix+=incx;
                     iy+=incy;
                  }
                  y[jy]+=temp*A[j]+alpha*sum;
                  A+=ldA;
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
                  y[j]+=temp*A[j];
                  for(i=j+1;i<n;i++)
                  {
                     y[i]+=temp*A[i];
                     sum+=A[i]*x[i];
                  }
                  y[j]+=alpha*sum;
                  A+=ldA;
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
                  y[jy]+=temp*A[j];
                  ix=jx;
                  iy=jy;
                  for(i=j+1;i<n;i++)
                  {
                     ix+=incx;
                     iy+=incy;
                     y[iy]+=temp*A[i];
                     sum+=A[i]*x[ix];
                  }
                  y[jy]+=alpha*sum;
                  jx+=incx;
                  jy+=incy;
                  A+=ldA;
               }
            }
         }      
      }
      return 0;
   }
   
   /// @brief Performs complex symmetric matrix-vector multiplication.
   ///
   /// For complex symmetric matrix A, complex vectors x and y, and complex scalars alpha and beta,
   ///
   ///         y := alpha*A*x + beta*y
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetric matrix A
   /// is to be referenced:  
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular,
   ///        if uplo = 'L' or 'l' then A is lower triangular.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param alpha Complex scalar.
   /// @param A Pointer to complex n-by-n symmetric matrix A.  
   /// Only the upper or lower triangular part of A is referenced, depending on the value of uplo above.
   /// @param ldA Column length of the matrix A.  ldA>=n
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of the vector x.
   /// @param beta Complex scalar.
   /// @param y Pointer to complex vector y.
   /// @param incy Increment of the vector y.
   /// @ingroup MATV

   template <typename real_t>
   int SYMV(char uplo, int_t n, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *x, int_t incx, complex<real_t> beta, complex<real_t> *y, int_t incy)
   {
      return SYMV<complex<real_t> >(uplo,n,alpha,A,ldA,x,incx,beta,y,incy);
   }
   
}
#endif

