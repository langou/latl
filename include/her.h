//
//  her.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/23/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _her_h
#define _her_h

/// @file her.h Performs complex vector outer product.


#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs a vector outer product of a complex vector with itself.
   /// 
   /// For a complex Hermitian matrix A, complex vector x and real scalar alpha,
   ///
   ///        A := alpha * x * x' + A
   ///
   /// is computed.  The result is Hermitian and is stored as either upper or lower triangular in A.
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
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of the vector x.  x!=0
   /// @param A Pointer to complex Hermitian n-by-n matrix A.  If uplo = 'U' or 'u' then only the upper
   /// triangular part of A is referenced and the lower part is not referenced.  If uplo = 'L' or 'l' then
   /// only the lower triangular part of A is referenced and the upper part is not referenced.
   /// @param ldA Column length of matrix A.  ldA>=n.
   /// @ingroup BLAS

   template <typename real_t>
   int HER(char uplo, int_t n, real_t alpha, complex<real_t> *x, int_t incx, complex<real_t> *A, int_t ldA)
   {
      using std::conj;
      using std::real;
      using std::toupper;
      const real_t zero(0.0);
      complex<real_t> t;
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
               t=alpha*conj(x[j]);
               for(i=0;i<j;i++)
                  A[i]+=x[i]*t;
               A[j]=real(A[j])+real(x[j]*t);
               A+=ldA;
            }
         }
         else
         {
            for(j=0;j<n;j++)
            {
               t=alpha*conj(x[j]);
               A[j]=real(A[j])+real(t*x[j]);
               for(i=j+1;i<n;i++)
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
               t=alpha*conj(x[jx]);
               ix=kx;
               for(i=0;i<j;i++)
               {
                  A[i]+=x[ix]*t;
                  ix+=incx;
               }
               A[j]=real(A[j])+real(x[jx]*t);
               A+=ldA;
               jx+=incx;
            }
         }
         else
         {
            jx=kx;
            for(j=0;j<n;j++)
            {
               t=alpha*conj(x[jx]);
               A[j]=real(A[j])+real(t*x[jx]);
               ix=jx;
               for(i=j+1;i<n;i++)
               {
                  ix+=incx;
                  A[i]+=x[ix]*t;
               }
               jx+=incx;
               A+=ldA;
            }
         }
      }
      return 0;
   }

#ifdef __latl_cblas
#include <cblas.h>

   template <> int HER<float>(char uplo, int_t n, float alpha, complex<float> *x, int_t incx, complex<float> *A, int_t ldA)
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

      cblas_cher(CblasColMajor,Uplo,n,alpha,x,incx,A,ldA);

      return 0;
   }
   
   template <> int HER<double>(char uplo, int_t n, double alpha, complex<double> *x, int_t incx, complex<double> *A, int_t ldA)
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

      cblas_zher(CblasColMajor,Uplo,n,alpha,x,incx,A,ldA);

      return 0;
   }
#endif
}
#endif
