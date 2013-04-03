//
//  geru.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/23/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _geru_h
#define _geru_h

/// @file geru.h Performs complex vector outer product.


#include "latl.h"

namespace LATL
{
   /// @brief Performs a vector outer product of two complex vectors.
   /// 
   /// For a complex matrix A, complex vectors x and y, and complex scalar alpha,
   ///
   ///        A := alpha * x * y.' + A
   ///
   /// is computed. 
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param m Specifies the number of rows of the matrix A.  m>=0
   /// @param n Specifies the number of coumns of the matrix A.  n>=0
   /// @param alpha Complex scalar.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of the vector x.  x!=0
   /// @param y Pointer to complex vector y.
   /// @param incy Increment of the vector y.  y!=0
   /// @param A Pointer to complex m-by-n matrix A.
   /// @param ldA Column length of matrix A.  ldA>=m.
   /// @ingroup BLAS

   template <typename real_t>
   int GERU(int_t m, int_t n, complex<real_t> alpha, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy, complex<real_t> *A, int_t ldA)
   {
      using std::conj;

      const complex<real_t> zero(0.0,0.0);
      int_t i,j,kx,jy,ix;
     
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(incx==0)
         return -5;
      else if(incy==0)
         return -7;
      else if(ldA<m)
         return -9;
      else if((m==0)||(n==0)||(alpha==zero))
         return 0;
      
      jy=(incy>0)?0:(1-n)*incy;
      kx=(incx>0)?0:(1-m)*incx;
      
      if(incx==1)
      {
         for(j=0;j<n;j++)
         {
            for(i=0;i<m;i++)
               A[i]+=x[i]*alpha*y[jy];
            jy+=incy;
            A+=ldA;
         }
      }
      else
      {
         for(j=0;j<n;j++)
         {
            ix=kx;
            for(i=0;i<m;i++)
            {
               A[i]+=x[ix]*alpha*y[jy];
               ix+=incx;
            }
            A+=ldA;
            jy+=incy;
         }
      }
      return 0;
   }
#ifdef __latl_cblas
#include <cblas.h>

   template <> int GERU<float>(int_t m, int_t n, complex<float> alpha, complex<float> *x, int_t incx, complex<float> *y, int_t incy, complex<float> *A, int_t ldA)
   {
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(incx==0)
         return -5;
      else if(incy==0)
         return -7;
      else if(ldA<m)
         return -9;

      cblas_cgeru(CblasColMajor,m,n,&alpha,x,incx,y,incy,A,ldA);

      return 0;
   }

   template <> int GERU<double>(int_t m, int_t n, complex<double> alpha, complex<double> *x, int_t incx, complex<double> *y, int_t incy, complex<double> *A, int_t ldA)
   {
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(incx==0)
         return -5;
      else if(incy==0)
         return -7;
      else if(ldA<m)
         return -9;

      cblas_cgeru(CblasColMajor,m,n,&alpha,x,incx,y,incy,A,ldA);

      return 0;
   }
   
#endif
}
#endif
