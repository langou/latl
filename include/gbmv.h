//
//  gbmv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/12/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _gbmv_h
#define _gbmv_h

/// @file gbmv.h Performs banded matrix vector multiplication.


#include <algorithm>
#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs banded real matrix vector multiplication.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies wheather the transpose of A is to be used or not:
   ///
   ///        if trans = 'N' or 'n' then y := alpha * A  * x  + beta * y
   ///        if trans = 'T' or 't' then y := alpha * A' * x  + beta * y
   ///        if trans = 'C' or 'c' then y := alpha * A' * x  + beta * y
   ///
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param kl Number of subdiagonals of the matrix A.  kl >= 0
   /// @param ku Number of superdiagonals of the matrix A.  ku >= 0
   /// @param alpha Real scalar.
   /// @param A Banded m-by-n real matrix.  The matrix bands are stored as rows in the matrix A as follows.
   /// The diagonal is stored as row ku, starting in column 0.  The first superdiagonal is stored starting in 
   /// column 1 of row ku-1, and the second superdiagonal in column 2 of row ku-2, and so on.  The first subdiagonal
   /// is stored starting in column 0 of row ku+1, and the second subdiagonal in column 0 of row ku+2, and so on.
   /// As an example, consider the following banded 5-by-5 matrix with ku=2 and kl=1.  On the left is the matrix
   /// in standard form, and on the right the same matrix in banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h d a . . ]        [ . . a b c ]
   ///        [ p i e b . ]        [ . d e f g ]
   ///        [ . q j f c ]        [ h i j k l ]
   ///        [ . . r k g ]        [ p q r s . ]
   ///        [ . . . s l ]
   ///
   /// @param lda Column length of matrix A.  lda >= kl+ku+1
   /// @param x Real vector of length n if trans='N' or 'n', or length m otherwise.
   /// @param incx Increment of vector x.
   /// @param beta Real scalar.
   /// @param y Real vector of length m if trans='N' or 'n', or length n otherwise.
   /// @param incy Increment of vector y.
   /// @ingroup MATV
   
   template <typename real_t>
   int GBMV(char trans, int_t m, int_t n, int_t kl, int_t ku, real_t alpha, real_t *A, int_t lda, real_t *x, int_t incx, real_t beta, real_t *y, int_t incy)
   {
      using std::min;
      using std::max;
      using std::toupper;
      const real_t one(1.0);
      const real_t zero(0.0);
      const int_t izero=0;
      int_t kx,ky,lenx,leny,i,ix,iy,j,jx,jy,k,i0,im;
      real_t temp;
      
      trans=toupper(trans);
      if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -1;
      else if(m<0)
         return -2;
      else if(n<0)
         return -3;
      else if(kl<0)
         return -4;
      else if(ku<0)
         return -5;
      else if(lda<kl+ku+1)
         return -8;
      else if(incx==0)
         return -10;
      else if(incy==0)
         return -13;
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
               for(j=0;j<n;j++)
               {
                  temp=alpha*x[jx];
                  k=ku-j;
                  i0=max(izero,j-ku);
                  im=min(m-1,j+kl);
                  for(i=i0;i<=im;i++)
                     y[i]+=temp*A[k+i];
                  A+=lda;
                  jx+=incx;
               }
            else
               for(j=0;j<n;j++)
               {
                  temp=alpha*x[jx];
                  iy=ky;
                  k=ku-j;
                  i0=max(izero,j-ku);
                  im=min(m-1,j+kl);
                  for(i=i0;i<=im;i++)
                  {
                     y[iy]+=temp*A[k+i];
                     iy+=incy;
                  }
                  A+=lda;
                  jx+=incx;
                  if(j>=ku)
                     ky+=incy;
               }
         }
         else
         {
            jy=ky;
            if(incx==1)
               for(j=0;j<n;j++)
               {
                  temp=zero;
                  k=ku-j;
                  i0=max(izero,j-ku);
                  im=min(m-1,j+kl);
                  for(i=i0;i<=im;i++)
                     temp+=A[k+i]*x[i];
                  y[jy]+=alpha*temp;
                  jy+=incy;
                  A+=lda;
               }
            else
               for(j=0;j<n;j++)
               {
                  temp=zero;
                  k=ku-j;
                  i0=max(izero,j-ku);
                  im=min(m-1,j+kl);
                  ix=kx;
                  for(i=i0;i<=im;i++)
                  {
                     temp+=A[k+i]*x[ix];
                     ix+=incx;
                  }
                  y[jy]+=alpha*temp;
                  jy+=incy;
                  A+=lda;
                  if(j>=ku)
                     kx+=incx;
               }
         }      
      }
      return 0;
   }
   
   /// @brief Performs banded complex matrix vector multiplication.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies wheather the transpose or conjugate transpose of A, or neither is used: 
   ///
   ///        if trans = 'N' or 'n' then y := alpha * A   * x + beta * y
   ///        if trans = 'T' or 't' then y := alpha * A.' * x + beta * y
   ///        if trans = 'C' or 'c' then y := alpha * A'  * x + beta * y
   ///
   /// @param m Number of rows of the matrix A.  m >= 0
   /// @param n Number of columns of the matrix A.  n >= 0
   /// @param kl Number of subdiagonals of the matrix A.  kl >= 0
   /// @param ku Number of superdiagonals of the matrix A.  ku >= 0
   /// @param alpha Complex scalar.
   /// @param A Banded m-by-n complex matrix.  The matrix bands are stored as rows in the matrix A as follows.
   /// The diagonal is stored as row ku, starting in column 0.  The first superdiagonal is stored starting in 
   /// column 1 of row ku-1, and the second superdiagonal in column 2 of row ku-2, and so on.  The first subdiagonal
   /// is stored starting in column 0 of row ku+1, and the second subdiagonal in column 0 of row ku+2, and so on.
   /// As an example, consider the following banded 5-by-5 matrix with ku=2 and kl=1.  On the left is the matrix
   /// in standard form, and on the right the same matrix in banded form.
   ///
   ///          (standard)           (banded)
   ///        [ h d a . . ]        [ . . a b c ]
   ///        [ p i e b . ]        [ . d e f g ]
   ///        [ . q j f c ]        [ h i j k l ]
   ///        [ . . r k g ]        [ p q r s . ]
   ///        [ . . . s l ]
   ///
   /// @param lda Column length of matrix A.  lda >= kl+ku+1
   /// @param x Complex vector of length n if trans='N' or 'n', or length m otherwise.
   /// @param incx Increment of vector x.
   /// @param beta Complex scalar.
   /// @param y Complex vector of length m if trans='N' or 'n', or length n otherwise.
   /// @param incy Increment of vector y.
   /// @ingroup MATV

   template <typename real_t>
   int GBMV(char trans, int_t m, int_t n, int_t kl, int_t ku, complex<real_t> alpha, complex<real_t> *A, int_t lda, complex<real_t> *x, int_t incx, complex<real_t> beta, complex<real_t> *y, int_t incy)
   {
      using std::conj;
      using std::min;
      using std::max;
      using std::toupper;
      const complex<real_t> one(1.0,0.0);
      const complex<real_t> zero(0.0,0.0);
      const int_t izero=0;
      int_t kx,ky,lenx,leny,i,ix,iy,j,jx,jy,k,i0,im;
      complex<real_t> temp;
      
      trans=toupper(trans);
      if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -1;
      else if(m<0)
         return -2;
      else if(n<0)
         return -3;
      else if(kl<0)
         return -4;
      else if(ku<0)
         return -5;
      else if(lda<kl+ku+1)
         return -8;
      else if(incx==0)
         return -10;
      else if(incy==0)
         return -13;
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
                     k=ku-j;
                     i0=max(izero,j-ku);
                     im=min(m-1,j+kl);
                     for(i=i0;i<=im;i++)
                        y[i]+=temp*A[k+i];
                     A+=lda;
                     jx+=incx;
                  }
               else
                  for(j=0;j<n;j++)
                  {
                     temp=alpha*x[jx];
                     iy=ky;
                     k=ku-j;
                     i0=max(izero,j-ku);
                     im=min(m-1,j+kl);
                     for(i=i0;i<=im;i++)
                     {
                        y[iy]+=temp*A[k+i];
                        iy+=incy;
                     }
                     A+=lda;
                     jx+=incx;
                     if(j>=ku)
                        ky+=incy;
                  }
               break;
               
            case 'T':
               jy=ky;
               if(incx==1)
                  for(j=0;j<n;j++)
                  {
                     temp=zero;
                     k=ku-j;
                     i0=max(izero,j-ku);
                     im=min(m-1,j+kl);
                     for(i=i0;i<=im;i++)
                        temp+=A[k+i]*x[i];
                     y[jy]+=alpha*temp;
                     jy+=incy;
                     A+=lda;
                  }
               else
                  for(j=0;j<n;j++)
                  {
                     temp=zero;
                     k=ku-j;
                     i0=max(izero,j-ku);
                     im=min(m-1,j+kl);
                     ix=kx;
                     for(i=i0;i<=im;i++)
                     {
                        temp+=A[k+i]*x[ix];
                        ix+=incx;
                     }
                     y[jy]+=alpha*temp;
                     jy+=incy;
                     A+=lda;
                     if(j>=ku)
                        kx+=incx;
                  }
               break;
               
            case 'C':
               jy=ky;
               if(incx==1)
                  for(j=0;j<n;j++)
                  {
                     temp=zero;
                     k=ku-j;
                     i0=max(izero,j-ku);
                     im=min(m-1,j+kl);
                     for(i=i0;i<=im;i++)
                        temp+=conj(A[k+i])*x[i];
                     y[jy]+=alpha*temp;
                     jy+=incy;
                     A+=lda;
                  }
               else
                  for(j=0;j<n;j++)
                  {
                     temp=zero;
                     k=ku-j;
                     i0=max(izero,j-ku);
                     im=min(m-1,j+kl);
                     ix=kx;
                     for(i=i0;i<=im;i++)
                     {
                        temp+=conj(A[k+i])*x[ix];
                        ix+=incx;
                     }
                     y[jy]+=alpha*temp;
                     jy+=incy;
                     A+=lda;
                     if(j>=ku)
                        kx+=incx;
                  }
               break;
         }      
      }
      return 0;
   }

}
#endif

