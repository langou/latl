//
//  trsm.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/5/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _trsm_h
#define _trsm_h

/// @file trsm.h Solves triangular matrix equations.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Solves a real triangular matrix equation.
   /// 
   /// For a real upper or lower triangular matrix A, real rectangular matrix B and real scalar alpha,
   ///
   ///        A*X=alpha*B  or  A'*X=alpha*B  or  X*A=alpha*B  or  X*A'=alpha*B
   /// is solved for the matrix X.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param side Specifies whether the matrix A appears on the left or right side as follows:
   ///
   ///        if side = 'L' or 'l' then A*X=alpha*B  or  A'*X=alpha*B is solved
   ///        if side = 'R' or 'r' then X*A=alpha*B  or  X*A'=alpha*B is solved
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param trans Specifies whether the transpose or conjugate transpose of A is to be used:
   ///
   ///        if trans = 'N' or 'n' then A*X=alpha*B or X*A=alpha*B is solved
   ///        if trans = 'T' or 't' then A'*X=alpha*B or X*A'=alpha*B is solved
   ///        if trans = 'C' or 'c' then A'*X=alpha*B or X*A'=alpha*B is solved
   /// @param diag specifies whether or not A is unit triangular as follows:
   ///
   ///        if diag = 'U' or 'u' then A is assumed to be unit triangular
   ///        if diag = 'N' or 'n' then A is not assumed to be unit triangular.
   /// If A is unit triangular, the diagonal elements are assumed to be unity and are not referenced.
   /// @param m Specifies the number of rows of the matrix B.  m>=0
   /// @param n Specifies the number of columns of the matrix B.  n>=0
   /// @param alpha Real scalar.
   /// @param A Pointer to real triangular matrix A.  The order of A is n if side = 'L' or 'l';
   /// the order of A is m is side = 'R' or 'r'.
   /// If uplo = 'U' or 'u, A is upper triangular and the lower triangular part is not referenced.  
   /// If uplo = 'L' or 'l A is lower triangular and the upper triangular part is not referenced.
   /// @param ldA Column length of the matrix A.  If side='L' or 'l' then ldA>=n; if side='R' or 'r' then ldA>=m.
   /// @param B Pointer to real m-by-n matrix B.  On entry, B contains the right hand side matrix B.  
   /// On exit, B is overwritten with the solution matrix X.
   /// @param ldB Column length of the matrix B.  ldB>=m.
   /// @ingroup BLAS

   template <typename real_t>
   int TRSM(char side, char uplo, char trans, char diag, int_t m, int_t n, real_t alpha, real_t *A, int_t ldA, real_t *B, int_t ldB)
   {
      using std::toupper;

      const real_t zero(0.0);
      const real_t one(1.0);
      int_t i,j,k;
      real_t *a,*b,*bt;
      real_t t;
      bool nounit;

      side=toupper(side);
      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);

      nounit=(diag=='N')?1:0;
      
      if((side!='L')&&(side!='R'))
         return -1;
      else if((uplo!='U')&&(uplo!='L'))
         return -2;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -3;
      else if((diag!='U')&&(diag!='N'))
         return -4;
      else if(m<0)
         return -5;
      else if(n<0)
         return -6;
      else if(ldA<((side=='L')?m:n))
         return -9;
      else if(ldB<m)
         return -11;
      else if((m==0)||(n==0))
         return 0;
      
      if(alpha==zero)
      {
         b=B;
         for(j=0;j<n;j++)
         {
            for(i=0;i<m;i++)
               b[i]=zero;
            b+=ldB;
         }
      }
      else
      {
         if(side=='L')
         {
            if(trans=='N')
            {
               if(uplo=='U')
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     a=A+m*ldA;
                     for(k=m-1;k>=0;k--)
                     {
                        a-=ldA;
                        if(nounit)
                           b[k]=b[k]/a[k];
                        for(i=0;i<k;i++)
                           b[i]-=b[k]*a[i];
                     }
                     b+=ldB;
                  }
               }
               else
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     a=A;
                     for(k=0;k<m;k++)
                     {
                        if(nounit)
                           b[k]=b[k]/a[k];
                        for(i=k+1;i<m;i++)
                           b[i]-=b[k]*a[i];
                        a+=ldA;
                     }
                     b+=ldB;
                  }
               }
            }
            else
            {
               if(uplo=='U')
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     a=A;
                     for(i=0;i<m;i++)
                     {
                        t=alpha*b[i];
                        for(k=0;k<i;k++)
                           t-=a[k]*b[k];
                        if(nounit)
                           t=t/a[i];
                        b[i]=t;
                        a+=ldA;
                     }
                     b+=ldB;
                  }
               }
               else
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     a=A+m*ldA;
                     for(i=m-1;i>=0;i--)
                     {
                        a-=ldA;
                        t=alpha*b[i];
                        for(k=i+1;k<m;k++)
                           t-=a[k]*b[k];
                        if(nounit)
                           t=t/a[i];
                        b[i]=t;
                     }
                     b+=ldB;
                  }
               }
            }
         }
         else
         {
            if(trans=='N')
            {
               if(uplo=='U')
               {
                  b=B;
                  a=A;
                  for(j=0;j<n;j++)
                  {
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     bt=B;
                     for(k=0;k<j;k++)
                     {
                        for(i=0;i<m;i++)
                           b[i]-=a[k]*bt[i];
                        bt+=ldB;
                     }
                     if(nounit)
                     {
                        t=one/a[j];
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                     a+=ldA;
                     b+=ldB;
                  }
               }
               else
               {
                  b=B+n*ldB;
                  a=A+n*ldA;
                  for(j=n-1;j>=0;j--)
                  {
                     a-=ldA;
                     b-=ldB;
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     bt=B+(j+1)*ldB;
                     for(k=j+1;k<n;k++)
                     {
                        for(i=0;i<m;i++)
                           b[i]-=a[k]*bt[i];
                        bt+=ldB;
                     }
                     if(nounit)
                     {
                        t=one/a[j];
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                  }
               }
            }
            else
            {
               if(uplo=='U')
               {
                  b=B+n*ldB;
                  a=A+n*ldA;
                  for(k=n-1;k>=0;k--)
                  {
                     a-=ldA;
                     b-=ldB;
                     if(nounit)
                     {
                        t=one/a[k];
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                     bt=B;
                     for(j=0;j<k;j++)
                     {
                        t=a[j];
                        for(i=0;i<m;i++)
                           bt[i]-=t*b[i];
                        bt+=ldB;
                     }
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                  }
               }
               else
               {
                  a=A;
                  b=B;
                  for(k=0;k<n;k++)
                  {
                     if(nounit)
                     {
                        t=one/a[k];
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                     bt=B+(k+1)*ldB;
                     for(j=k+1;j<n;j++)
                     {
                        t=a[j];
                        for(i=0;i<m;i++)
                           bt[i]-=t*b[i];
                        bt+=ldB;
                     }
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     a+=ldA;
                     b+=ldB;
                  }
               }
            }
         }
      }
      return 0;
   }

   /// @brief Solves a complex triangular matrix equation.
   /// 
   /// For a complex upper or lower triangular matrix A, complex rectangular matrix B and complex scalar alpha,
   ///
   ///        A*X=alpha*B or A'*X=alpha*B or A.'*X=alpha*B or X*A=alpha*B or X*A'=alpha*B or X*A.'=alpha*B
   /// is solved for the matrix X.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param side Specifies whether the matrix A appears on the left or right side as follows:
   ///
   ///        if side = 'L' or 'l' then A*X=alpha*B  or  A'*X=alpha*B is solved
   ///        if side = 'R' or 'r' then X*A=alpha*B  or  X*A'=alpha*B is solved
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param trans Specifies whether the transpose of A is to be used or not:
   ///
   ///        if trans = 'N' or 'n' then A*X=alpha*B or X*A=alpha*B is solved
   ///        if trans = 'T' or 't' then A.'*X=alpha*B or X*A.'=alpha*B is solved
   ///        if trans = 'C' or 'c' then A'*X=alpha*B or X*A'=alpha*B is solved
   /// @param diag specifies whether or not A is unit triangular as follows:
   ///
   ///        if diag = 'U' or 'u' then A is assumed to be unit triangular
   ///        if diag = 'N' or 'n' then A is not assumed to be unit triangular.
   /// If A is unit triangular, the diagonal elements are assumed to be unity and are not referenced.
   /// @param m Specifies the number of rows of the matrix B.  m>=0
   /// @param n Specifies the number of columns of the matrix B.  n>=0
   /// @param alpha Complex scalar.
   /// @param A Pointer to complex triangular matrix A.  The order of A is n if side = 'L' or 'l';
   /// the order of A is m is side = 'R' or 'r'.
   /// If uplo = 'U' or 'u, A is upper triangular and the lower triangular part is not referenced.  
   /// If uplo = 'L' or 'l A is lower triangular and the upper triangular part is not referenced.
   /// @param ldA Column length of the matrix A.  If side='L' or 'l' then ldA>=n; if side='R' or 'r' then ldA>=m.
   /// @param B Pointer to complex m-by-n matrix B.  On entry, B contains the right hand side matrix B.  
   /// On exit, B is overwritten with the solution matrix X.
   /// @param ldB Column length of the matrix B.  ldB>=m.
   /// @ingroup BLAS

   template <typename real_t>
   int TRSM(char side, char uplo, char trans, char diag, int_t m, int_t n, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *B, int_t ldB)
   {
      
      using std::conj;
      using std::toupper;
      const complex<real_t> zero(0.0,0.0);
      const complex<real_t> one(1.0,0.0);
      int_t i,j,k;
      complex<real_t> *a,*b,*bt;
      complex<real_t> t;
      bool nounit;
      
      side=toupper(side);
      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);

      nounit=(diag=='N')?1:0;
      
      if((side!='L')&&(side!='R'))
         return -1;
      else if((uplo!='U')&&(uplo!='L'))
         return -2;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -3;
      else if((diag!='U')&&(diag!='N'))
         return -4;
      else if(m<0)
         return -5;
      else if(n<0)
         return -6;
      else if(ldA<((side=='L')?m:n))
         return -9;
      else if(ldB<m)
         return -11;
      else if((m==0)||(n==0))
         return 0;
      
      if(alpha==zero)
      {
         b=B;
         for(j=0;j<n;j++)
         {
            for(i=0;i<m;i++)
               b[i]=zero;
            b+=ldB;
         }
      }
      else
      {
         if(side=='L')
         {
            if(trans=='N')
            {
               if(uplo=='U')
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     a=A+m*ldA;
                     for(k=m-1;k>=0;k--)
                     {
                        a-=ldA;
                        if(nounit)
                           b[k]=b[k]/a[k];
                        for(i=0;i<k;i++)
                           b[i]-=b[k]*a[i];
                     }
                     b+=ldB;
                  }
               }
               else
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     a=A;
                     for(k=0;k<m;k++)
                     {
                        if(nounit)
                           b[k]=b[k]/a[k];
                        for(i=k+1;i<m;i++)
                           b[i]-=b[k]*a[i];
                        a+=ldA;
                     }
                     b+=ldB;
                  }
               }
            }
            else if(trans=='T')
            {
               if(uplo=='U')
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     a=A;
                     for(i=0;i<m;i++)
                     {
                        t=alpha*b[i];
                        for(k=0;k<i;k++)
                           t-=a[k]*b[k];
                        if(nounit)
                           t=t/a[i];
                        b[i]=t;
                        a+=ldA;
                     }
                     b+=ldB;
                  }
               }
               else
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     a=A+m*ldA;
                     for(i=m-1;i>=0;i--)
                     {
                        a-=ldA;
                        t=alpha*b[i];
                        for(k=i+1;k<m;k++)
                           t-=a[k]*b[k];
                        if(nounit)
                           t=t/a[i];
                        b[i]=t;
                     }
                     b+=ldB;
                  }
               }
            }
            else if(trans=='C')
            {
               if(uplo=='U')
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     a=A;
                     for(i=0;i<m;i++)
                     {
                        t=alpha*b[i];
                        for(k=0;k<i;k++)
                           t-=conj(a[k])*b[k];
                        if(nounit)
                           t=t/conj(a[i]);
                        b[i]=t;
                        a+=ldA;
                     }
                     b+=ldB;
                  }
               }
               else
               {
                  b=B;
                  for(j=0;j<n;j++)
                  {
                     a=A+m*ldA;
                     for(i=m-1;i>=0;i--)
                     {
                        a-=ldA;
                        t=alpha*b[i];
                        for(k=i+1;k<m;k++)
                           t-=conj(a[k])*b[k];
                        if(nounit)
                           t=t/conj(a[i]);
                        b[i]=t;
                     }
                     b+=ldB;
                  }
               }
            }
         }
         else
         {
            if(trans=='N')
            {
               if(uplo=='U')
               {
                  b=B;
                  a=A;
                  for(j=0;j<n;j++)
                  {
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     bt=B;
                     for(k=0;k<j;k++)
                     {
                        for(i=0;i<m;i++)
                           b[i]-=a[k]*bt[i];
                        bt+=ldB;
                     }
                     if(nounit)
                     {
                        t=one/a[j];
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                     a+=ldA;
                     b+=ldB;
                  }
               }
               else
               {
                  b=B+n*ldB;
                  a=A+n*ldA;
                  for(j=n-1;j>=0;j--)
                  {
                     a-=ldA;
                     b-=ldB;
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     bt=B+(j+1)*ldB;
                     for(k=j+1;k<n;k++)
                     {
                        for(i=0;i<m;i++)
                           b[i]-=a[k]*bt[i];
                        bt+=ldB;
                     }
                     if(nounit)
                     {
                        t=one/a[j];
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                  }
               }
            }
            else if(trans=='T')
            {
               if(uplo=='U')
               {
                  b=B+n*ldB;
                  a=A+n*ldA;
                  for(k=n-1;k>=0;k--)
                  {
                     a-=ldA;
                     b-=ldB;
                     if(nounit)
                     {
                        t=one/a[k];
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                     bt=B;
                     for(j=0;j<k;j++)
                     {
                        t=a[j];
                        for(i=0;i<m;i++)
                           bt[i]-=t*b[i];
                        bt+=ldB;
                     }
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                  }
               }
               else
               {
                  a=A;
                  b=B;
                  for(k=0;k<n;k++)
                  {
                     if(nounit)
                     {
                        t=one/a[k];
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                     bt=B+(k+1)*ldB;
                     for(j=k+1;j<n;j++)
                     {
                        t=a[j];
                        for(i=0;i<m;i++)
                           bt[i]-=t*b[i];
                        bt+=ldB;
                     }
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     a+=ldA;
                     b+=ldB;
                  }
               }
            }
            else if(trans=='C')
            {
               if(uplo=='U')
               {
                  b=B+n*ldB;
                  a=A+n*ldA;
                  for(k=n-1;k>=0;k--)
                  {
                     a-=ldA;
                     b-=ldB;
                     if(nounit)
                     {
                        t=one/conj(a[k]);
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                     bt=B;
                     for(j=0;j<k;j++)
                     {
                        t=conj(a[j]);
                        for(i=0;i<m;i++)
                           bt[i]-=t*b[i];
                        bt+=ldB;
                     }
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                  }
               }
               else
               {
                  a=A;
                  b=B;
                  for(k=0;k<n;k++)
                  {
                     if(nounit)
                     {
                        t=one/conj(a[k]);
                        for(i=0;i<m;i++)
                           b[i]*=t;
                     }
                     bt=B+(k+1)*ldB;
                     for(j=k+1;j<n;j++)
                     {
                        t=conj(a[j]);
                        for(i=0;i<m;i++)
                           bt[i]-=t*b[i];
                        bt+=ldB;
                     }
                     for(i=0;i<m;i++)
                        b[i]*=alpha;
                     a+=ldA;
                     b+=ldB;
                  }
               }
            }
         }
      }
      return 0;
   }

#ifdef __latl_cblas
#include <cblas.h>

   template <> int TRSM<float>(char side, char uplo, char trans, char diag, int_t m, int_t n, float alpha, float *A, int_t ldA, float *B, int_t ldB)
   {
      using std::toupper;
      side=toupper(side);
      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);
      if((side!='L')&&(side!='R'))
         return -1;
      else if((uplo!='U')&&(uplo!='L'))
         return -2;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -3;
      else if((diag!='U')&&(diag!='N'))
         return -4;
      else if(m<0)
         return -5;
      else if(n<0)
         return -6;
      else if(ldA<((side=='L')?m:n))
         return -9;
      else if(ldB<m)
         return -11;
      else if((m==0)||(n==0))
         return 0;

      const CBLAS_SIDE Side=(side=='L')?CblasLeft:CblasRight;
      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);
      const CBLAS_DIAG Diag=(diag=='N')?CblasNonUnit:CblasUnit;

      cblas_strsm(CblasColMajor,Side,Uplo,Trans,Diag,m,n,alpha,A,ldA,B,ldB);

      return 0;
   }

   template <> int TRSM<double>(char side, char uplo, char trans, char diag, int_t m, int_t n, double alpha, double *A, int_t ldA, double *B, int_t ldB)
   {
      using std::toupper;
      side=toupper(side);
      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);
      if((side!='L')&&(side!='R'))
         return -1;
      else if((uplo!='U')&&(uplo!='L'))
         return -2;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -3;
      else if((diag!='U')&&(diag!='N'))
         return -4;
      else if(m<0)
         return -5;
      else if(n<0)
         return -6;
      else if(ldA<((side=='L')?m:n))
         return -9;
      else if(ldB<m)
         return -11;
      else if((m==0)||(n==0))
         return 0;

      const CBLAS_SIDE Side=(side=='L')?CblasLeft:CblasRight;
      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);
      const CBLAS_DIAG Diag=(diag=='N')?CblasNonUnit:CblasUnit;

      cblas_dtrsm(CblasColMajor,Side,Uplo,Trans,Diag,m,n,alpha,A,ldA,B,ldB);

      return 0;
   }
   template <> int TRSM<float>(char side, char uplo, char trans, char diag, int_t m, int_t n, complex<float> alpha, complex<float> *A, int_t ldA, complex<float> *B, int_t ldB)
   {
      using std::toupper;
      side=toupper(side);
      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);
      if((side!='L')&&(side!='R'))
         return -1;
      else if((uplo!='U')&&(uplo!='L'))
         return -2;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -3;
      else if((diag!='U')&&(diag!='N'))
         return -4;
      else if(m<0)
         return -5;
      else if(n<0)
         return -6;
      else if(ldA<((side=='L')?m:n))
         return -9;
      else if(ldB<m)
         return -11;
      else if((m==0)||(n==0))
         return 0;

      const CBLAS_SIDE Side=(side=='L')?CblasLeft:CblasRight;
      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);
      const CBLAS_DIAG Diag=(diag=='N')?CblasNonUnit:CblasUnit;

      cblas_ctrsm(CblasColMajor,Side,Uplo,Trans,Diag,m,n,&alpha,A,ldA,B,ldB);

      return 0;
   }

   template <> int TRSM<double>(char side, char uplo, char trans, char diag, int_t m, int_t n, complex<double> alpha, complex<double> *A, int_t ldA, complex<double> *B, int_t ldB)
   {
      using std::toupper;
      side=toupper(side);
      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);
      if((side!='L')&&(side!='R'))
         return -1;
      else if((uplo!='U')&&(uplo!='L'))
         return -2;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -3;
      else if((diag!='U')&&(diag!='N'))
         return -4;
      else if(m<0)
         return -5;
      else if(n<0)
         return -6;
      else if(ldA<((side=='L')?m:n))
         return -9;
      else if(ldB<m)
         return -11;
      else if((m==0)||(n==0))
         return 0;

      const CBLAS_SIDE Side=(side=='L')?CblasLeft:CblasRight;
      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);
      const CBLAS_DIAG Diag=(diag=='N')?CblasNonUnit:CblasUnit;

      cblas_ztrsm(CblasColMajor,Side,Uplo,Trans,Diag,m,n,&alpha,A,ldA,B,ldB);
      
      return 0;
   }
#endif
}
#endif

