//
//  gemv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/12/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _gemv_h
#define _gemv_h

/// @file gemv.h Computes general matrix-vector products.

#include "latl.h"

namespace latl
{
   /// @brief Computes general real matrix-vector products.
   ///
   /// For a real matrix A and real scalars alpha, beta
   ///
   ///        y := alpha * op(A) * x + beta * y
   ///
   /// is computed, where op(A) = A or A'.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies whether the transpose of A is used or not:
   ///
   ///        if trans = 'N' or 'n' then op(A) = A
   ///        if trans = 'T' or 't' then op(A) = A'
   ///        if trans = 'C' or 'c' then op(A) = A'
   ///
   /// @param m The number of rows of the matrix A.  m>=0
   /// @param n The number of columns of the matrix A.  n>=0
   /// @param alpha Real scalar.
   /// @param A Pointer to real matrix m-by-n matrix A.  
   /// @param ldA Column length of the matrix A.  ldA>=m.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of the vector x.  
   /// @param beta Real scalar.
   /// @param y Pointer to real vector y.  
   /// @param incy Increment of the vector y.
   /// @ingroup MATV

   template <typename real_t>
   int gemv(char trans, int_t m, int_t n, real_t alpha, real_t *A, int_t ldA, real_t *x, int_t incx, real_t beta, real_t *y, int_t incy)
   {
      const real_t one(1.0);
      const real_t zero(0.0);
      int_t kx,ky,lenx,leny,i,ix,iy,j,jx,jy;
      real_t temp;
      
      if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -1;
      else if(m<0)
         return -2;
      else if(n<0)
         return -3;
      else if(ldA<m)
         return -6;
      else if(incx==0)
         return -8;
      else if(incy==0)
         return -11;
      else if((m==0)||(n==0))
         return 0;
      
      lenx=(trans=='N')?n:m;
      leny=(trans=='N')?m:n;
      kx=(incx>0)?0:(1-lenx)*incx;
      ky=(incy>0)?0:(1-leny)*incy;
      
      if(beta==zero)
      {
         if(incy==1)
         {
            for(i=0;i<leny;i++)
               y[i]=zero;
         }
         else
         {
            iy=ky;
            for(i=0;i<leny;i++)
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
            for(i=0;i<leny;i++)
               y[i]*=beta;
         }
         else
         {
            iy=ky;
            for(i=0;i<leny;i++)
            {
               y[iy]*=beta;
               iy+=incy;
            }
         }
      }
      
      if(alpha!=zero)
      {
         if(trans=='N')
         {
            jx=kx;
            if(incy==1)
            {
               for(j=0;j<n;j++)
               {
                  temp=alpha*x[jx];
                  for(i=0;i<m;i++)
                     y[i]+=temp*A[i];
                  A+=ldA;
                  jx+=incx;
               }
            }
            else
            {
               for(j=0;j<n;j++)
               {
                  temp=alpha*x[jx];
                  iy=ky;
                  for(i=0;i<m;i++)
                  {
                     y[iy]+=temp*A[i];
                     iy+=incy;
                  }
                  A+=ldA;
                  jx+=incx;
               }
            }
         }
         else
         {
            jy=ky;
            if(incx==1)
            {
               for(j=0;j<n;j++)
               {
                  temp=zero;
                  for(i=0;i<m;i++)
                     temp+=A[i]*x[i];
                  y[jy]+=alpha*temp;
                  jy+=incy;
                  A+=ldA;
               }
            }
            else
            {
               for(j=0;j<n;j++)
               {
                  temp=zero;
                  ix=kx;
                  for(i=0;i<m;i++)
                  {
                     temp+=A[i]*x[ix];
                     ix+=incx;
                  }
                  y[jy]+=alpha*temp;
                  jy+=incy;
                  A+=ldA;
               }
            }
         }      
      }
      return 0;
   }
   
   /// @brief Computes general complex matrix-vector products.
   ///
   /// For a complex matrix A and complex scalars alpha, beta
   ///
   ///        y := alpha * op(A) * x + beta * y
   ///
   /// is computed, where op(A) = A or A'.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies whether the transpose or conjugate transpose of A, or neither is used:
   ///
   ///        if trans = 'N' or 'n' then op(A) = A
   ///        if trans = 'T' or 't' then op(A) = A.'
   ///        if trans = 'C' or 'c' then op(A) = A'
   ///
   /// @param m The number of rows of the matrix A.  m>=0
   /// @param n The number of columns of the matrix A.  n>=0
   /// @param alpha Complex scalar.
   /// @param A Pointer to complex matrix m-by-n matrix A.  
   /// @param ldA Column length of the matrix A.  ldA>=m.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of the vector x.  
   /// @param beta Complex scalar.
   /// @param y Pointer to complex vector y.  
   /// @param incy Increment of the vector y.
   /// @ingroup MATV

   template <typename real_t>
   int gemv(char trans, int_t m, int_t n, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *x, int_t incx, complex<real_t> beta, complex<real_t> *y, int_t incy)
   {
      using std::conj;
      const complex<real_t> one(1.0,0.0);
      const complex<real_t> zero(0.0,0.0);
      int_t kx,ky,lenx,leny,i,ix,iy,j,jx,jy;
      complex<real_t> temp;
      
      if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -1;
      else if(m<0)
         return -2;
      else if(n<0)
         return -3;
      else if(ldA<m)
         return -6;
      else if(incx==0)
         return -8;
      else if(incy==0)
         return -11;
      else if((m==0)||(n==0))
         return 0;
      
      lenx=(trans=='N')?n:m;
      leny=(trans=='N')?m:n;
      kx=(incx>0)?0:(1-lenx)*incx;
      ky=(incy>0)?0:(1-leny)*incy;
      
      if(beta==zero)
      {
         if(incy==1)
         {
            for(i=0;i<leny;i++)
               y[i]=zero;
         }
         else
         {
            iy=ky;
            for(i=0;i<leny;i++)
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
            for(i=0;i<leny;i++)
               y[i]*=beta;
         }
         else
         {
            iy=ky;
            for(i=0;i<leny;i++)
            {
               y[iy]*=beta;
               iy+=incy;
            }
         }
      }
      
      if(alpha!=zero)
      {
         switch(trans)
         {
            case 'N':
               jx=kx;
               if(incy==1)
                  for(j=0;j<n;j++)
                  {
                     temp=alpha*x[jx];
                     for(i=0;i<m;i++)
                        y[i]+=temp*A[i];
                     A+=ldA;
                     jx+=incx;
                  }
               else
                  for(j=0;j<n;j++)
                  {
                     temp=alpha*x[jx];
                     iy=ky;
                     for(i=0;i<m;i++)
                     {
                        y[iy]+=temp*A[i];
                        iy+=incy;
                     }
                     A+=ldA;
                     jx+=incx;
                  }
               break;
               
            case 'T':
               jy=ky;
               if(incx==1)
                  for(j=0;j<n;j++)
                  {
                     temp=zero;
                     for(i=0;i<m;i++)
                        temp+=A[i]*x[i];
                     y[jy]+=alpha*temp;
                     jy+=incy;
                     A+=ldA;
                  }
               else
                  for(j=0;j<n;j++)
                  {
                     temp=zero;
                     ix=kx;
                     for(i=0;i<m;i++)
                     {
                        temp+=A[i]*x[ix];
                        ix+=incx;
                     }
                     y[jy]+=alpha*temp;
                     jy+=incy;
                     A+=ldA;
                  }
               
               break;
               
            case 'C':
               jy=ky;
               if(incx==1)
                  for(j=0;j<n;j++)
                  {
                     temp=zero;
                     for(i=0;i<m;i++)
                        temp+=conj(A[i])*x[i];
                     y[jy]+=alpha*temp;
                     jy+=incy;
                     A+=ldA;
                  }
               else
                  for(j=0;j<n;j++)
                  {
                     temp=zero;
                     ix=kx;
                     for(i=0;i<m;i++)
                     {
                        temp+=conj(A[i])*x[ix];
                        ix+=incx;
                     }
                     y[jy]+=alpha*temp;
                     jy+=incy;
                     A+=ldA;
                  }
               
               break;
         }
      }
      return 0;
   }
   
}
#endif
