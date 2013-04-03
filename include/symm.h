//
//  symm.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/2/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _symm_h
#define _symm_h

/// @file symm.h Performs general symmetric matrix-matrix multiplication.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief  Performs general real symmetric matrix-matrix multiplication.
   ///
   /// For real matrices B and C, symmetric matrix A, and real scalars alpha and beta,
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
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetric matrix A
   /// is to be referenced:  
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular,
   ///        if uplo = 'L' or 'l' then A is lower triangular.
   /// @param m Specifies the number of rows of the matrices B and C.  m>=0
   /// @param n Specifies the number of columns of the matrices B and C.  n>=0
   /// @param alpha Real scalar.
   /// @param A Pointer to real symmetric matrix A.  If side = 'L' or 'l', then A is m-by-m;
   /// if side = 'R' or 'r', then A is n-by-n.
   /// @param ldA Column length of the matrix A.  If side = 'L' or 'l', ldA>=m.  If side = 'R' or 'r', ldA>=n.
   /// @param B Pointer to real m-by-n matrix B.
   /// @param ldB Column length of the matrix B.  ldB>=m.
   /// @param beta Real scalar.
   /// @param C Pointer to real m-by-n matrix C.
   /// @param ldC Column length of the matrix C.  ldC>=m.
   /// @ingroup BLAS

   template <typename real_t>
   int SYMM(char side, char uplo, int_t m, int_t n, real_t alpha, real_t *A, int_t ldA, real_t *B, int_t ldB, real_t beta, real_t *C, int_t ldC)
   {
      using std::toupper;

      const real_t zero(0.0);
      const real_t one(1.0);
      int_t i,j,k;
      real_t *a,*b,*c,*at,*bt;
      real_t s,t;

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
                     s+=b[k]*a[k];
                  }
                  c[i]=beta*c[i]+t*a[i]+alpha*s;
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
                     s+=b[k]*a[k];
                  }
                  c[i]=beta*c[i]+t*a[i]+alpha*s;
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
               t=alpha*a[j];
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
                  t=alpha*at[j];
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
               t=alpha*a[j];
               for(i=0;i<m;i++)
                  c[i]=c[i]*beta+t*b[i];
               bt=B;
               at=A;
               for(k=0;k<j;k++)
               {
                  t=alpha*at[j];
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
   
   /// @brief  Performs general complex symmetric matrix-matrix multiplication.
   ///
   /// For complex matrices B and C, symmetric matrix A, and complex scalars alpha and beta,
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
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetric matrix A
   /// is to be referenced:  
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular,
   ///        if uplo = 'L' or 'l' then A is lower triangular.
   /// @param m Specifies the number of rows of the matrices B and C.  m>=0
   /// @param n Specifies the number of columns of the matrices B and C.  n>=0
   /// @param alpha Complex scalar.
   /// @param A Pointer to complex symmetric matrix A.  If side = 'L' or 'l', then A is m-by-m;
   /// if side = 'R' or 'r', then A is n-by-n.
   /// @param ldA Column length of the matrix A.  If side = 'L' or 'l', ldA>=m.  If side = 'R' or 'r', ldA>=n.
   /// @param B Pointer to complex m-by-n matrix B.
   /// @param ldB Column length of the matrix B.  ldB>=m.
   /// @param beta Complex scalar.
   /// @param C Pointer to complex m-by-n matrix C.
   /// @param ldC Column length of the matrix C.  ldC>=m.
   /// @ingroup BLAS

   template <typename real_t>
   int SYMM(char side, char uplo, int_t m, int_t n, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *B, int_t ldB, complex<real_t> beta, complex<real_t> *C, int_t ldC)
   {
      return SYMM< complex<real_t> >(side,uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC);
   }

#ifdef __latl_cblas
#include <cblas.h>
   template <> int SYMM<float>(char side, char uplo, int_t m, int_t n, float alpha, float *A, int_t ldA, float *B, int_t ldB, float beta, float *C, int_t ldC)
   {
      using std::toupper;
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

      const CBLAS_SIDE Side=(side=='L')?CblasLeft:CblasRight;
      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;

      cblas_ssymm(CblasColMajor,Side,Uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC);

      return 0;

   }
   
   template <> int SYMM<double>(char side, char uplo, int_t m, int_t n, double alpha, double *A, int_t ldA, double *B, int_t ldB, double beta, double *C, int_t ldC)
   {
      using std::toupper;
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

      const CBLAS_SIDE Side=(side=='L')?CblasLeft:CblasRight;
      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;

      cblas_dsymm(CblasColMajor,Side,Uplo,m,n,alpha,A,ldA,B,ldB,beta,C,ldC);

      return 0;
      
   }
   template <> int SYMM<float>(char side, char uplo, int_t m, int_t n, complex<float> alpha, complex<float> *A, int_t ldA, complex<float> *B, int_t ldB, complex<float> beta, complex<float> *C, int_t ldC)
   {
      using std::toupper;
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

      const CBLAS_SIDE Side=(side=='L')?CblasLeft:CblasRight;
      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;

      cblas_csymm(CblasColMajor,Side,Uplo,m,n,&alpha,A,ldA,B,ldB,&beta,C,ldC);

      return 0;
      
   }
   template <> int SYMM<double>(char side, char uplo, int_t m, int_t n, complex<double> alpha, complex<double> *A, int_t ldA, complex<double> *B, int_t ldB, complex<double> beta, complex<double> *C, int_t ldC)
   {
      using std::toupper;
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

      const CBLAS_SIDE Side=(side=='L')?CblasLeft:CblasRight;
      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;

      cblas_zsymm(CblasColMajor,Side,Uplo,m,n,&alpha,A,ldA,B,ldB,&beta,C,ldC);

      return 0;
      
   }
#endif
}
#endif

