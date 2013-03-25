//
//  trmm.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/4/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _trmm_h
#define _trmm_h

/// @file trmm.h Performs multiplication of a triangular matrix with a rectangular matrix.

#include <cctype>
#include "latl.h"

namespace LATL
{   
   /// @brief Performs multiplication of a real triangular matrix with a real rectangular matrix.
   /// 
   /// For a real upper or lower triangular matrix A, real rectangular matrix B and real scalar alpha,
   ///
   ///        B := alpha*A*B  or  B := alpha*A'*B  or  B := alpha*B*A  or  B := alpha*B*A'
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param side Specifies whether the matrix A appears on the left or right side as follows:
   ///
   ///        if side = 'L' or 'l' then B := alpha*A*B or alpha*A'*B
   ///        if side = 'R' or 'r' then C := alpha*B*A or alpha*B*A'
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param trans Specifies wheather the transpose of A is to be used or not:
   ///
   ///        if trans = 'N' or 'n' then B := alpha*A*B or alpha*B*A
   ///        if trans = 'T' or 't' then B := alpha*A'*B or alpha*B*A'
   ///        if trans = 'C' or 'c' then B := alpha*A'*B or alpha*B*A'
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
   /// @param ldA Column length of the matrix A.  If side='L' or 'l' then ldA>=n; if side='R' or 'r' then lda>=m.
   /// @param B Pointer to real m-by-n matrix B.
   /// @param ldB Column length of the matrix B.  ldB>=m.
   /// @ingroup MATM

   template <typename real_t>
   int TRMM(char side, char uplo, char trans, char diag, int_t m, int_t n, real_t alpha, real_t *A, int_t ldA, real_t *B, int_t ldB)
   {
      using std::toupper;
      
      const real_t zero(0.0);
      int_t i,j,k;
      real_t *a,*b,*bt;
      real_t t;
      
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

      bool nounit=(diag=='N')?1:0;

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
      else if(side=='L')
      {
         if(trans=='N')
         {
            if(uplo=='U')
            {
               b=B;
               for(j=0;j<n;j++)
               {
                  a=A;
                  for(k=0;k<m;k++)
                  {
                     t=alpha*b[k];
                     for(i=0;i<k;i++)
                        b[i]+=t*a[i];
                     if(nounit)
                        t*=a[k];
                     b[k]=t;
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
                  for(k=m-1;k>=0;k--)
                  {
                     a-=ldA;
                     t=alpha*b[k];
                     b[k]=t;
                     if(nounit)
                        b[k]*=a[k];
                     for(i=k+1;i<m;i++)
                        b[i]+=t*a[i];
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
                  a=A+m*ldA;
                  for(i=m-1;i>=0;i--)
                  {
                     a-=ldA;
                     t=b[i];
                     if(nounit)
                        t*=a[i];
                     for(k=0;k<i;k++)
                        t+=a[k]*b[k];
                     b[i]=alpha*t;
                  }
                  b+=ldB;
               }
            }
            else
            {
               b=B;
               for(j=0;j<n;j++)
               {
                  a=A;
                  for(i=0;i<m;i++)
                  {
                     t=b[i];
                     if(nounit)
                        t*=a[i];
                     for(k=i+1;k<m;k++)
                        t+=a[k]*b[k];
                     b[i]=alpha*t;
                     a+=ldA;
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
               a=A+n*ldA;
               b=B+n*ldB;
               for(j=n-1;j>=0;j--)
               {
                  a-=ldA;
                  b-=ldB;
                  if(nounit)
                     t=alpha*a[j];
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     b[i]*=t;
                  bt=B;
                  for(k=0;k<j;k++)
                  {
                     t=alpha*a[k];
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     bt+=ldB;
                  }
               }
            }
            else
            {
               a=A;
               b=B;
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     t=alpha*a[j];
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     b[i]*=t;
                  bt=b+ldB;
                  for(k=j+1;k<n;k++)
                  {
                     t=alpha*a[k];
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     bt+=ldB;
                  }
                  a+=ldA;
                  b+=ldB;
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               a=A;
               bt=B;
               for(k=0;k<n;k++)
               {
                  b=B;
                  for(j=0;j<k;j++)
                  {
                     t=alpha*a[j];
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     b+=ldB;
                  }
                  if(nounit)
                     t=alpha*a[k];
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     b[i]*=t;
                  a+=ldA;
                  bt+=ldB;
               }
            }
            else
            {
               a=A+n*ldA;
               bt=B+n*ldB;
               for(k=n-1;k>=0;k--)
               {
                  a-=ldA;
                  bt-=ldB;
                  b=B+(k+1)*ldB;
                  for(j=k+1;j<n;j++)
                  {
                     t=alpha*a[j];
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     b+=ldB;
                  }
                  if(nounit)
                     t=alpha*a[k];
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     bt[i]*=t;
               }
            }
         }
      }
      return 0;
   }

   /// @brief Performs multiplication of a complex triangular matrix with a complex rectangular matrix.
   /// 
   /// For a complex upper or lower triangular matrix A, complex rectangular matrix B and complex scalar alpha,
   ///
   ///        B:=alpha*A*B or B:=alpha*A'*B or B:=alpha*A.'*B or B:=alpha*B*A or B:=alpha*B*A' or B:=alpha*B*A.'
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param side Specifies whether the matrix A appears on the left or right side as follows:
   ///
   ///        if side = 'L' or 'l' then B := alpha*A*B or alpha*A'*B or alpha*A.'*B
   ///        if side = 'R' or 'r' then C := alpha*B*A or alpha*B*A' or alpha*B*A.'
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param trans Specifies wheather the transpose or conjugate transpose is applied to A as follows:
   ///
   ///        if trans = 'N' or 'n' then B := alpha*A*B or alpha*B*A
   ///        if trans = 'T' or 't' then B := alpha*A.'*B or alpha*B*A.'
   ///        if trans = 'C' or 'c' then B := alpha*A'*B or alpha*B*A'
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
   /// @param ldA Column length of the matrix A.  If side='L' or 'l' then ldA>=n; if side='R' or 'r' then lda>=m.
   /// @param B Pointer to complex m-by-n matrix B.
   /// @param ldB Column length of the matrix B.  ldB>=m.
   /// @ingroup MATM

   template <typename real_t>
   int TRMM(char side, char uplo, char trans, char diag, int_t m, int_t n, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *B, int_t ldB)
   {
      using std::conj;
      using std::toupper;
      
      const complex<real_t> zero(0.0,0.0);
      int_t i,j,k;
      complex<real_t> *a,*b,*bt;
      complex<real_t> t;
      
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

      bool nounit=(diag=='N')?1:0;

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
      else if(side=='L')
      {
         if(trans=='N')
         {
            if(uplo=='U')
            {
               b=B;
               for(j=0;j<n;j++)
               {
                  a=A;
                  for(k=0;k<m;k++)
                  {
                     t=alpha*b[k];
                     for(i=0;i<k;i++)
                        b[i]+=t*a[i];
                     if(nounit)
                        t*=a[k];
                     b[k]=t;
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
                  for(k=m-1;k>=0;k--)
                  {
                     a-=ldA;
                     t=alpha*b[k];
                     b[k]=t;
                     if(nounit)
                        b[k]*=a[k];
                     for(i=k+1;i<m;i++)
                        b[i]+=t*a[i];
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
                  a=A+m*ldA;
                  for(i=m-1;i>=0;i--)
                  {
                     a-=ldA;
                     t=b[i];
                     if(nounit)
                        t*=a[i];
                     for(k=0;k<i;k++)
                        t+=a[k]*b[k];
                     b[i]=alpha*t;
                  }
                  b+=ldB;
               }
            }
            else
            {
               b=B;
               for(j=0;j<n;j++)
               {
                  a=A;
                  for(i=0;i<m;i++)
                  {
                     t=b[i];
                     if(nounit)
                        t*=a[i];
                     for(k=i+1;k<m;k++)
                        t+=a[k]*b[k];
                     b[i]=alpha*t;
                     a+=ldA;
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
                  a=A+m*ldA;
                  for(i=m-1;i>=0;i--)
                  {
                     a-=ldA;
                     t=b[i];
                     if(nounit)
                        t*=conj(a[i]);
                     for(k=0;k<i;k++)
                        t+=conj(a[k])*b[k];
                     b[i]=alpha*t;
                  }
                  b+=ldB;
               }
            }
            else
            {
               b=B;
               for(j=0;j<n;j++)
               {
                  a=A;
                  for(i=0;i<m;i++)
                  {
                     t=b[i];
                     if(nounit)
                        t*=conj(a[i]);
                     for(k=i+1;k<m;k++)
                        t+=conj(a[k])*b[k];
                     b[i]=alpha*t;
                     a+=ldA;
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
               a=A+n*ldA;
               b=B+n*ldB;
               for(j=n-1;j>=0;j--)
               {
                  a-=ldA;
                  b-=ldB;
                  if(nounit)
                     t=alpha*a[j];
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     b[i]*=t;
                  bt=B;
                  for(k=0;k<j;k++)
                  {
                     t=alpha*a[k];
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     bt+=ldB;
                  }
               }
            }
            else
            {
               a=A;
               b=B;
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     t=alpha*a[j];
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     b[i]*=t;
                  bt=b+ldB;
                  for(k=j+1;k<n;k++)
                  {
                     t=alpha*a[k];
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     bt+=ldB;
                  }
                  a+=ldA;
                  b+=ldB;
               }
            }
         }
         else if(trans=='T')
         {
            if(uplo=='U')
            {
               a=A;
               bt=B;
               for(k=0;k<n;k++)
               {
                  b=B;
                  for(j=0;j<k;j++)
                  {
                     t=alpha*a[j];
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     b+=ldB;
                  }
                  if(nounit)
                     t=alpha*a[k];
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     b[i]*=t;
                  a+=ldA;
                  bt+=ldB;
               }
            }
            else
            {
               a=A+n*ldA;
               bt=B+n*ldB;
               for(k=n-1;k>=0;k--)
               {
                  a-=ldA;
                  bt-=ldB;
                  b=B+(k+1)*ldB;
                  for(j=k+1;j<n;j++)
                  {
                     t=alpha*a[j];
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     b+=ldB;
                  }
                  if(nounit)
                     t=alpha*a[k];
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     bt[i]*=t;
               }
            }
         }
         else if(trans=='C')
         {
            if(uplo=='U')
            {
               a=A;
               bt=B;
               for(k=0;k<n;k++)
               {
                  b=B;
                  for(j=0;j<k;j++)
                  {
                     t=alpha*conj(a[j]);
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     b+=ldB;
                  }
                  if(nounit)
                     t=alpha*conj(a[k]);
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     b[i]*=t;
                  a+=ldA;
                  bt+=ldB;
               }
            }
            else
            {
               a=A+n*ldA;
               bt=B+n*ldB;
               for(k=n-1;k>=0;k--)
               {
                  a-=ldA;
                  bt-=ldB;
                  b=B+(k+1)*ldB;
                  for(j=k+1;j<n;j++)
                  {
                     t=alpha*conj(a[j]);
                     for(i=0;i<m;i++)
                        b[i]+=t*bt[i];
                     b+=ldB;
                  }
                  if(nounit)
                     t=alpha*conj(a[k]);
                  else
                     t=alpha;
                  for(i=0;i<m;i++)
                     bt[i]*=t;
               }
            }
         }
      }
      return 0;
   }

#ifdef __latl_cblas
#include <cblas.h>

   template <> int TRMM<float>(char side, char uplo, char trans, char diag, int_t m, int_t n, float alpha, float *A, int_t ldA, float *B, int_t ldB)
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
      const CBLAS_DIAG Diag=(diag=='N')?CblasNonUnit:CblasUnit;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:CblasTrans;

      cblas_strmm(CblasColMajor,Side,Uplo,Trans,Diag,m,n,alpha,A,ldA,B,ldB);
      
      return 0;
   }

   template <> int TRMM<double>(char side, char uplo, char trans, char diag, int_t m, int_t n, double alpha, double *A, int_t ldA, double *B, int_t ldB)
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
      const CBLAS_DIAG Diag=(diag=='N')?CblasNonUnit:CblasUnit;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:CblasTrans;

      cblas_dtrmm(CblasColMajor,Side,Uplo,Trans,Diag,m,n,alpha,A,ldA,B,ldB);

      return 0;
   }
   template <> int TRMM<float>(char side, char uplo, char trans, char diag, int_t m, int_t n, complex<float> alpha, complex<float> *A, int_t ldA, complex<float> *B, int_t ldB)
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
      const CBLAS_DIAG Diag=(diag=='N')?CblasNonUnit:CblasUnit;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);

      cblas_ctrmm(CblasColMajor,Side,Uplo,Trans,Diag,m,n,&alpha,A,ldA,B,ldB);

      return 0;
   }
   template <> int TRMM<double>(char side, char uplo, char trans, char diag, int_t m, int_t n, complex<double> alpha, complex<double> *A, int_t ldA, complex<double> *B, int_t ldB)
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
      const CBLAS_DIAG Diag=(diag=='N')?CblasNonUnit:CblasUnit;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);

      cblas_ztrmm(CblasColMajor,Side,Uplo,Trans,Diag,m,n,&alpha,A,ldA,B,ldB);

      return 0;
   }
#endif

}
#endif
