//
//  spr.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/23/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _spr_h
#define _spr_h

/// @file spr.h Performs vector outer products using packed format.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs a real vector outer product using packed storage.
   /// 
   /// For real vector x and real scalar alpha
   ///
   ///        A := alpha*x*x'+A
   /// is computed, where the real symmetric matrix A uses packed storage format.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the packed matrix A is upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular,
   ///        if uplo = 'L' or 'l' then A is lower triangular.
   /// @param n Specifies the order of the symmetric matrix A.  n>=0
   /// @param alpha Real scalar.
   /// @param x Pointer the real vector x.
   /// @param incx Increment of vector x. incx!=0
   /// @param A Pointer to packed real symmetric n-by-n matrix A.  The packed storage format stores a triangular matrix
   /// using variable length columns so that no space is wasted.  If A is upper triangular, the first column consists only
   /// of the diagonal element, which has length 1.  The second column has length 2, the third length 3, and so on until the
   /// nth column, which has length n.  Similarly, if A is lower triangular, the first column has length n, the second length n-1,
   /// and so on until the nth column, which has length 1.  The entire n-by-n triangular matrix is stored using n(n+1)/2 elements.
   /// @ingroup BLAS

   template <typename real_t>
   int SPR(char uplo, int_t n, real_t alpha, real_t *x, int_t incx, real_t *A)
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
               A+=j+1;
            }
         }
         else
         {
            for(j=0;j<n;j++)
            {
               t=alpha*x[j];
               for(i=j;i<n;i++)
                  A[i-j]+=x[i]*t;
               A+=n-j;
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
               A+=j+1;
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
                  A[i-j]+=x[ix]*t;
                  ix+=incx;
               }
               jx+=incx;
               A+=n-j;
            }
         }
      }
      return 0;
   }
   
   /// @brief Performs a complex vector outer product using packed storage.
   /// 
   /// For complex vector x and complex scalar alpha
   ///
   ///        A := alpha*x*x.'+A
   /// is computed, where the complex symmetric matrix A uses packed storage format.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the packed matrix A is upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular,
   ///        if uplo = 'L' or 'l' then A is lower triangular.
   /// @param n Specifies the order of the symmetric matrix A.  n>=0
   /// @param alpha Complex scalar.
   /// @param x Pointer the complex vector x.
   /// @param incx Increment of vector x. incx!=0
   /// @param A Pointer to packed complex symmetric n-by-n matrix A.  The packed storage format stores a triangular matrix
   /// using variable length columns so that no space is wasted.  If A is upper triangular, the first column consists only
   /// of the diagonal element, which has length 1.  The second column has length 2, the third length 3, and so on until the
   /// nth column, which has length n.  Similarly, if A is lower triangular, the first column has length n, the second length n-1,
   /// and so on until the nth column, which has length 1.  The entire n-by-n triangular matrix is stored using n(n+1)/2 elements.
   /// @ingroup BLAS

   template <typename real_t>
   int SPR(char uplo, int_t n, complex<real_t> alpha, complex<real_t> *x, int_t incx, complex<real_t> *A)
   {
      return SPR<complex<real_t> >(uplo,n,alpha,x,incx,A);
   }
}
#endif
