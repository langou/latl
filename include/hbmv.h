//
//  hbmv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/12/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _hbmv_h
#define _hbmv_h

/// @file hbmv.h Performs banded Hermitian complex matrix-vector multiplication.


#include <algorithm>
#include <cctype>
#include "latl.h"

namespace LATL
{

   /// @brief Performs banded complex Hermitian matrix-vector multiplication.
   /// 
   /// For a banded complex Hermitian matrix A, complex vectors x and y, and real scalars alpha and beta,
   ///
   ///        y := alpha * A * x + beta * y
   ///
   /// is computed.  
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   ///
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param k Specifies the number of super-diagonals of the matrix A.  k >= 0
   /// @param alpha Complex scalar.
   /// @param A Pointer to banded complex Hermitian n-by-n matrix A.  The bands of A are stored as rows,
   /// while preserving the columns of A.  If uplo = 'U' or 'u' then only the upper triangular part of A 
   /// is referenced and the lower part is not referenced; the diagonal is stored on row k, the first super-diagonal
   /// in row k-1 starting in column 1, the second super-diagonal in row k-2 starting in column 2, and so on.  
   /// As an example, consider the following Hermitian matrix with n=5 and k=2.  On the left is the matrix in
   /// standard upper triangular form, and on the right is the same matrix in upper triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h d a . . ]        [ . . a b c ]
   ///        [ . i e b . ]        [ . d e f g ]
   ///        [ . . j f c ]        [ h i j k l ]
   ///        [ . . . k g ]        
   ///        [ . . . . l ]
   /// If uplo = 'L' or 'l' then only the lower triangular part of A is referenced and the upper part is not referenced;
   /// the diagonal is stored in row zero, the first sub-diagonal in row 1, and second sub-diagonal in row 2, and so on.  
   /// As an example, consider the following Hermitian matrix with n=5 and k=2.  On the left is the matrix in standard 
   /// lower triangular form, and on the right in the same matrix in lower triangular banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h . . . . ]        [ h i j k l ]
   ///        [ d i . . . ]        [ d e f g . ]
   ///        [ a e j . . ]        [ a b c . . ]
   ///        [ . b f k . ]        
   ///        [ . . c g l ]
   /// @param ldA Column length of matrix A.  ldA>=k+1.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of the vector x.  incx!=0
   /// @param beta Complex scalar.
   /// @param y Pointer to complex vector y.
   /// @param incy Increment of the vector y.  incy!=0
   /// @ingroup MATV
   
   template <typename real_t>
   int hbmv(char uplo, int_t n, int_t k, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *x, int_t incx, complex<real_t> beta, complex<real_t> *y, int_t incy)
   {
      using std::conj;
      using std::real;
      using std::max;
      using std::min;
      using std::toupper;
      
      const complex<real_t> one(1.0,0.0);
      const complex<real_t> zero(0.0,0.0);
      int_t kx,ky,i,ix,iy,j,jx,jy,i0,in;
      const int_t izero=0;

      complex<real_t> temp,sum;

      uplo=toupper(uplo);

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(k<0)
         return -3;
      else if(ldA<k+1)
         return -6;
      else if(incx==0)
         return -8;
      else if(incy==0)
         return -11;
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
                  i0=max(izero,j-k);
                  for(i=i0;i<j;i++)
                  {
                     y[i]+=temp*A[k-j+i];
                     sum+=conj(A[k-j+i])*x[i];
                  }
                  y[j]+=temp*real(A[k])+alpha*sum;
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
                  i0=max(izero,j-k);
                  for(i=i0;i<j;i++)
                  {
                     y[iy]+=temp*A[k-j+i];
                     sum+=conj(A[k-j+i])*x[ix];
                     ix+=incx;
                     iy+=incy;
                  }
                  y[jy]+=temp*real(A[k])+alpha*sum;
                  A+=ldA;
                  jx+=incx;
                  jy+=incy;
                  if(j>=k)
                  {
                     kx+=incx;
                     ky+=incy;
                  }
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
                  in=min(n-1,j+k);
                  for(i=j+1;i<=in;i++)
                  {
                     y[i]+=temp*A[i-j];
                     sum+=conj(A[i-j])*x[i];
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
                  ix=jx;
                  iy=jy;
                  in=min(n-1,j+k);
                  y[jy]+=temp*real(A[0]);
                  for(i=j+1;i<=in;i++)
                  {
                     ix+=incx;
                     iy+=incy;
                     y[iy]+=temp*A[i-j];
                     sum+=conj(A[i-j])*x[ix];
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
}
#endif
