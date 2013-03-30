//
//  ger.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/23/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _ger_h
#define _ger_h

/// @file ger.h Performs vector outer product.


#include "latl.h"

namespace LATL
{
   /// @brief Performs a vector outer product of two real vectors.
   /// 
   /// For a real matrix A, real vectors x and y, and real scalar alpha,
   ///
   ///        A := alpha * x * y' + A
   ///
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param m Specifies the number of rows of the matrix A.  m>=0
   /// @param n Specifies the number of coumns of the matrix A.  n>=0
   /// @param alpha Real scalar.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of the vector x.  x!=0
   /// @param y Pointer to real vector y.
   /// @param incy Increment of the vector y.  y!=0
   /// @param A Pointer to real m-by-n matrix A.
   /// @param ldA Column length of matrix A.  ldA>=m.
   /// @ingroup BLAS

   template <typename real_t>
   int GER(int_t m, int_t n, real_t alpha, real_t *x, int_t incx, real_t *y, int_t incy, real_t *A, int_t ldA)
   {
      const real_t zero(0.0);
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
   int GER(int_t m, int_t n, complex<real_t> alpha, complex<real_t> *x, int_t incx, complex<real_t> *y, int_t incy, complex<real_t> *A, int_t ldA)
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
}
#endif
