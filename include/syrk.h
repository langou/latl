//
//  syrk.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/4/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _syrk_h
#define _syrk_h

/// @file syrk.h Performs matrix-matrix multiplication operations resulting in symmetric matrices.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs multiplcation of a real matrix with its transpose.
   ///
   /// For real matrix A, real symmetric matrix C, and real scalars alpha and beta
   ///
   ///        C := alpha*A*A'+beta*C  or  C := alpha*A'*A+beta*C
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetric matrix C
   /// is to be referenced:  
   ///
   ///        if uplo = 'U' or 'u' then C is upper triangular,
   ///        if uplo = 'L' or 'l' then C is lower triangular.
   /// @param trans Specifies the operation to be perfomed as follows:
   ///
   ///        if trans = 'N' or 'n' then C := alpha*A*A'+beta*C
   ///        if trans = 'T' or 't' then C := alpha*A'*A+beta*C   
   ///        if trans = 'C' or 'c' then C := alpha*A'*A+beta*C
   /// @param n Specifies the order of the symmetric matrix C.  n>=0
   /// @param k Specifies the other dimension of the real matrix A (see below). k>=0
   /// @param alpha Real scalar.
   /// @param A Pointer to real matrix.
   ///
   ///        if trans = 'N' or 'n' then A is n-by-k
   ///        if trans = 'T' or 't' then A is k-by-n
   ///        if trans = 'C' or 'c' then A is k-by-n
   /// @param ldA Column length of the matrix A. If trans = 'N' or 'n' ldA>=n, otherwise ldA>=k.
   /// @param beta Real scalar.
   /// @param C Pointer to real symmetric n-by-n matrix C.  
   /// Only the upper or lower triangular part of C is referenced, depending on the value of uplo above.
   /// @param ldC Column length of the matrix C.  ldC>=n
   /// @ingroup BLAS

   template <typename real_t>
   int SYRK(char uplo, char trans, int_t n, int_t k, real_t alpha, real_t *A, int_t ldA, real_t beta, real_t *C, int_t ldC)
   {
      using std::toupper;

      const real_t zero(0.0);
      const real_t one(1.0);
      int_t i,j,l;
      real_t *a,*c,*at;
      real_t t;
      
      uplo=toupper(uplo);
      trans=toupper(trans);
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -2;
      else if(n<0)
         return -3;
      else if(k<0)
         return -4;
      else if(ldA<((trans=='N')?n:k))
         return -7;
      else if(ldC<n)
         return -10;
      else if((n==0)||(((alpha==zero)||(k==0))&&(beta==one)))
         return 0;
      
      if(alpha==zero)
      {
         if(uplo=='U')
         {
            c=C;
            for(j=0;j<n;j++)
            {
               for(i=0;i<=j;i++)
                  c[i]*=beta;
               c+=ldC;
            }
            
         }
         else
         {
            c=C;
            for(j=0;j<n;j++)
            {
               for(i=j;i<n;i++)
                  c[i]*=beta;
               c+=ldC;
            }
         }
      }
      else if(trans=='N')
      {
         if(uplo=='U')
         {
            c=C;
            for(j=0;j<n;j++)
            {
               for(i=0;i<=j;i++)
                  c[i]*=beta;
               a=A;
               for(l=0;l<k;l++)
               {
                  t=alpha*a[j];
                  for(i=0;i<=j;i++)
                     c[i]+=t*a[i];
                  a+=ldA;
               }
               c+=ldC;
            }
         }
         else
         {
            c=C;
            for(j=0;j<n;j++)
            {
               for(i=j;i<n;i++)
                  c[i]*=beta;
               a=A;
               for(l=0;l<k;l++)
               {
                  t=alpha*a[j];
                  for(i=j;i<n;i++)
                     c[i]+=t*a[i];
                  a+=ldA;
               }
               c+=ldC;
            }
         }
      }
      else 
      {
         if(uplo=='U')
         {
            c=C;
            at=A;
            for(j=0;j<n;j++)
            {
               a=A;
               for(i=0;i<=j;i++)
               {
                  t=zero;
                  for(l=0;l<k;l++)
                     t+=a[l]*at[l];
                  c[i]=alpha*t+beta*c[i];
                  a+=ldA;
               }
               at+=ldA;
               c+=ldC;
            }
         }
         else
         {
            at=A;
            c=C;
            for(j=0;j<n;j++)
            {
               a=A+j*ldA;
               for(i=j;i<n;i++)
               {
                  t=zero;
                  for(l=0;l<k;l++)
                     t+=a[l]*at[l];
                  c[i]=alpha*t+beta*c[i];
                  a+=ldA;
               }
               at+=ldA;
               c+=ldC;
            }
         }
      }
      return 0;
   }
   
   /// @brief Performs multiplcation of a complex matrix with its transpose.
   ///
   /// For complex matrix A, complex symmetric matrix C, and complex scalars alpha and beta
   ///
   ///        C := alpha*A*A.'+beta*C  or  C := alpha*A.'*A+beta*C
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetric matrix C
   /// is to be referenced:  
   ///
   ///        if uplo = 'U' or 'u' then C is upper triangular,
   ///        if uplo = 'L' or 'l' then C is lower triangular.
   /// @param trans Specifies the operation to be perfomed as follows:
   ///
   ///        if trans = 'N' or 'n' then C := alpha*A*A.'+beta*C
   ///        if trans = 'T' or 't' then C := alpha*A.'*A+beta*C
   /// @param n Specifies the order of the complex symmetric matrix C.  n>=0
   /// @param k Specifies the other dimension of the complex matrix A (see below). k>=0
   /// @param alpha Complex scalar.
   /// @param A Pointer to complex matrix.
   ///
   ///        if trans = 'N' or 'n' then A is n-by-k
   ///        if trans = 'T' or 't' then A is k-by-n
   /// @param ldA Column length of the matrix A. If trans = 'N' or 'n' ldA>=n, otherwise ldA>=k.
   /// @param beta Complex scalar.
   /// @param C Pointer to complex symmetric n-by-n matrix C.  
   /// Only the upper or lower triangular part of C is referenced, depending on the value of uplo above.
   /// @param ldC Column length of the matrix C.  ldC>=n
   /// @ingroup BLAS

   template <typename real_t>
   int SYRK(char uplo, char trans, int_t n, int_t k, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> beta, complex<real_t> *C, int_t ldC)
   {
      using std::toupper;

      const complex<real_t> zero(0.0,0.0);
      const complex<real_t> one(1.0,0.0);
      int_t i,j,l;
      complex<real_t> *a,*c,*at;
      complex<real_t> t;
      
      uplo=toupper(uplo);
      trans=toupper(trans);

      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T'))
         return -2;
      else if(n<0)
         return -3;
      else if(k<0)
         return -4;
      else if(ldA<((trans=='N')?n:k))
         return -7;
      else if(ldC<n)
         return -10;
      else if((n==0)||(((alpha==zero)||(k==0))&&(beta==one)))
         return 0;
      
      if(alpha==zero)
      {
         if(uplo=='U')
         {
            c=C;
            for(j=0;j<n;j++)
            {
               for(i=0;i<=j;i++)
                  c[i]*=beta;
               c+=ldC;
            }
            
         }
         else
         {
            c=C;
            for(j=0;j<n;j++)
            {
               for(i=j;i<n;i++)
                  c[i]*=beta;
               c+=ldC;
            }
         }
      }
      else if(trans=='N')
      {
         if(uplo=='U')
         {
            c=C;
            for(j=0;j<n;j++)
            {
               for(i=0;i<=j;i++)
                  c[i]*=beta;
               a=A;
               for(l=0;l<k;l++)
               {
                  t=alpha*a[j];
                  for(i=0;i<=j;i++)
                     c[i]+=t*a[i];
                  a+=ldA;
               }
               c+=ldC;
            }
         }
         else
         {
            c=C;
            for(j=0;j<n;j++)
            {
               for(i=j;i<n;i++)
                  c[i]*=beta;
               a=A;
               for(l=0;l<k;l++)
               {
                  t=alpha*a[j];
                  for(i=j;i<n;i++)
                     c[i]+=t*a[i];
                  a+=ldA;
               }
               c+=ldC;
            }
         }
      }
      else 
      {
         if(uplo=='U')
         {
            c=C;
            at=A;
            for(j=0;j<n;j++)
            {
               a=A;
               for(i=0;i<=j;i++)
               {
                  t=zero;
                  for(l=0;l<k;l++)
                     t+=a[l]*at[l];
                  c[i]=alpha*t+beta*c[i];
                  a+=ldA;
               }
               at+=ldA;
               c+=ldC;
            }
         }
         else
         {
            at=A;
            c=C;
            for(j=0;j<n;j++)
            {
               a=A+j*ldA;
               for(i=j;i<n;i++)
               {
                  t=zero;
                  for(l=0;l<k;l++)
                     t+=a[l]*at[l];
                  c[i]=alpha*t+beta*c[i];
                  a+=ldA;
               }
               at+=ldA;
               c+=ldC;
            }
         }
      }
      return 0;
   }

#ifdef __latl_cblas
#include <cblas.h>

   template <> int SYRK<float>(char uplo, char trans, int_t n, int_t k, float alpha, float *A, int_t ldA, float beta, float *C, int_t ldC)
   {
      using std::toupper;
      uplo=toupper(uplo);
      trans=toupper(trans);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -2;
      else if(n<0)
         return -3;
      else if(k<0)
         return -4;
      else if(ldA<((trans=='N')?n:k))
         return -7;
      else if(ldC<n)
         return -10;

      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);

      cblas_ssyrk(CblasColMajor,Uplo,Trans,n,k,alpha,A,ldA,beta,C,ldC);

      return 0;

   }

   template <> int SYRK<double>(char uplo, char trans, int_t n, int_t k, double alpha, double *A, int_t ldA, double beta, double *C, int_t ldC)
   {
      using std::toupper;
      uplo=toupper(uplo);
      trans=toupper(trans);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -2;
      else if(n<0)
         return -3;
      else if(k<0)
         return -4;
      else if(ldA<((trans=='N')?n:k))
         return -7;
      else if(ldC<n)
         return -10;

      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);

      cblas_dsyrk(CblasColMajor,Uplo,Trans,n,k,alpha,A,ldA,beta,C,ldC);

      return 0;
      
   }

   template <> int SYRK<float>(char uplo, char trans, int_t n, int_t k, complex<float> alpha, complex<float> *A, int_t ldA, complex<float> beta, complex<float> *C, int_t ldC)
   {
      using std::toupper;
      uplo=toupper(uplo);
      trans=toupper(trans);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T'))
         return -2;
      else if(n<0)
         return -3;
      else if(k<0)
         return -4;
      else if(ldA<((trans=='N')?n:k))
         return -7;
      else if(ldC<n)
         return -10;

      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);

      cblas_csyrk(CblasColMajor,Uplo,Trans,n,k,&alpha,A,ldA,&beta,C,ldC);

      return 0;

   }

   template <> int SYRK<double>(char uplo, char trans, int_t n, int_t k, complex<double> alpha, complex<double> *A, int_t ldA, complex<double> beta, complex<double> *C, int_t ldC)
   {
      using std::toupper;
      uplo=toupper(uplo);
      trans=toupper(trans);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T'))
         return -2;
      else if(n<0)
         return -3;
      else if(k<0)
         return -4;
      else if(ldA<((trans=='N')?n:k))
         return -7;
      else if(ldC<n)
         return -10;

      const CBLAS_UPLO Uplo=(uplo=='U')?CblasUpper:CblasLower;
      const CBLAS_TRANSPOSE Trans=(trans=='N')?CblasNoTrans:((trans=='T')?CblasTrans:CblasConjTrans);

      cblas_zsyrk(CblasColMajor,Uplo,Trans,n,k,&alpha,A,ldA,&beta,C,ldC);

      return 0;
      
   }

#endif
}

#endif

