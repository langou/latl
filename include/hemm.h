//
//  hemm.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/6/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _hemm_h
#define _hemm_h

/// @file hemm.h Performs general Hermitian complex matrix-matrix multiplication.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief  Performs general Hermitian complex matrix-matrix multiplication.
   ///
   /// For complex matrices B and C, Hermitian matrix A, and complex scalars alpha and beta,
   ///
   ///         C := alpha*A*B + beta*C  or   C := alpha*B*A + beta*C 
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param side Specifies whether the matrix A appears on the left or right side as follows:
   ///
   ///        if side = 'L' or 'l' then C := alpha*A*B + beta*C,
   ///        if side = 'R' or 'r' then C := alpha*B*A + beta*C.
   /// @param uplo Specifies whether the upper or lower triangular part of the Hermitian matrix A
   /// is to be referenced:  
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular,
   ///        if uplo = 'L' or 'l' then A is lower triangular.
   /// @param m Specifies the number of rows of the matrices B and C.  m>=0
   /// @param n Specifies the number of columns of the matrices B and C.  n>=0
   /// @param alpha Complex scalar.
   /// @param A Pointer to complex Hermitian matrix A.  If side = 'L' or 'l', then A is m-by-m;
   /// if side = 'R' or 'r', then A is n-by-n.
   /// @param ldA Column length of the matrix A.  If side = 'L' or 'l', ldA>=m.  If side = 'R' or 'r', ldA>=n.
   /// @param B Pointer to complex m-by-n matrix B.
   /// @param ldB Column length of the matrix B.  ldB>=m.
   /// @param beta Complex scalar.
   /// @param C Pointer to complex m-by-n matrix C.
   /// @param ldC Column length of the matrix C.  ldC>=m.
   /// @ingroup MATM

   template <typename real_t>
   int HEMM(char side, char uplo, int_t m, int_t n, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *B, int_t ldB, complex<real_t> beta, complex<real_t> *C, int_t ldC)
   {
      using std::conj;
      using std::real;
      using std::toupper;
      
      const complex<real_t> zero(0.0,0.0);
      const complex<real_t> one(1.0,0.0);
      int_t i,j,k;
      complex<real_t> *a,*b,*c,*at,*bt;
      complex<real_t> s,t;

      side=toupper(side);
      uplo=toupper(uplo);

      if((side!='L')&&(side!='R'))
         return -1;
      else if((uplo!='U')&&(uplo!='L'))
         return -2;
      else if(m<0)
         return -3;
      else if(n<0)
         return -4;
      else if(ldA<((side=='L')?m:n))
         return -7;
      else if(ldB<m)
         return -9;
      else if(ldC<m)
         return -12;
      else if((m==0)||(n==0)||((alpha==zero)&&(beta==one)))
         return 0;
      
      if(alpha==zero)
      {
         c=C;
         for(j=0;j<n;j++)
         {
            for(i=0;i<m;i++)
               c[i]*=beta;
            c+=ldC;
         }
      }
      else if(side=='L')
      {
         if(uplo=='U')
         {
            b=B;
            c=C;
            for(j=0;j<n;j++)
            {
               a=A;
               for(i=0;i<m;i++)
               {
                  t=alpha*b[i];
                  s=zero;
                  for(k=0;k<i;k++)
                  {
                     c[k]+=t*a[k];
                     s+=b[k]*conj(a[k]);
                  }
                  c[i]=beta*c[i]+t*real(a[i])+alpha*s;
                  a+=ldA;
               }
               b+=ldB;
               c+=ldC;
            }
         }
         else
         {
            b=B;
            c=C;
            for(j=0;j<n;j++)
            {
               a=A+m*ldA;
               for(i=m-1;i>=0;i--)
               {
                  a-=ldA;
                  t=alpha*b[i];
                  s=zero;
                  for(k=i+1;k<m;k++)
                  {
                     c[k]+=t*a[k];
                     s+=b[k]*conj(a[k]);
                  }
                  c[i]=beta*c[i]+t*real(a[i])+alpha*s;
               }
               b+=ldB;
               c+=ldC;
            }
         }
      }
      else 
      {
         if(uplo=='U')
         {
            a=A;
            c=C;
            b=B;
            for(j=0;j<n;j++)
            {
               t=alpha*real(a[j]);
               for(i=0;i<m;i++)
                  c[i]=c[i]*beta+t*b[i];
               bt=B;
               for(k=0;k<j;k++)
               {
                  t=alpha*a[k];
                  for(i=0;i<m;i++)
                     c[i]+=bt[i]*t;
                  bt+=ldB;
               }
               at=A+(j+1)*ldA;
               bt=B+(j+1)*ldB;
               for(k=j+1;k<n;k++)
               {
                  t=alpha*conj(at[j]);
                  for(i=0;i<m;i++)
                     c[i]+=t*bt[i];
                  at+=ldA;
                  bt+=ldB;
               }
               a+=ldA;
               b+=ldB;
               c+=ldC;
            }
         }
         else
         {
            a=A;
            c=C;
            b=B;
            for(j=0;j<n;j++)
            {
               t=alpha*real(a[j]);
               for(i=0;i<m;i++)
                  c[i]=c[i]*beta+t*b[i];
               bt=B;
               at=A;
               for(k=0;k<j;k++)
               {
                  t=alpha*conj(at[j]);
                  for(i=0;i<m;i++)
                     c[i]+=bt[i]*t;
                  at+=ldA;
                  bt+=ldB;
               }
               bt=B+(j+1)*ldB;
               for(k=j+1;k<n;k++)
               {
                  t=alpha*at[k];
                  for(i=0;i<m;i++)
                     c[i]+=t*bt[i];
                  bt+=ldB;
               }
               a+=ldA;
               b+=ldB;
               c+=ldC;
            }
         }
      }
      return 0;
   }
}
#endif
