//
//  syr2k.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/4/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _syr2k_h
#define _syr2k_h

/// @file syr2k.h Performs multiplication of two matrices resulting in a symmetric matrix.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs multiplcation of real matrices resulting in a symmetric matrix.
   ///
   /// For real matrices A and B, real symmetric matrix C, and real scalars alpha and beta
   ///
   ///        C := alpha*A*B'+alpha*B*A'+beta*C
   /// or
   ///
   ///        C := alpha*A'*B+alpha*B'*A+beta*C
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
   ///        if trans = 'N' or 'n' then C := alpha*A*B'+alpha*B*A'+beta*C
   ///        if trans = 'T' or 't' then C := alpha*A'*B+alpha*B'*A+beta*C
   ///        if trans = 'C' or 'c' then C := alpha*A'*B+alpha*B'*A+beta*C
   /// @param n Specifies the order of the complex Hermitian matrix C.  n>=0
   /// @param k Specifies the other dimension of the complex matrices A and B.
   ///
   ///        if trans = 'N' or 'n' then A and B are n-by-k
   ///        if trans = 'T' or 't' then A and B are k-by-n
   ///        if trans = 'C' or 'c' then A and B are k-by-n
   /// @param alpha Real scalar.
   /// @param A Pointer to real matrix.
   /// @param ldA Column length of the matrix A. If trans = 'N' or 'n' ldA>=n, otherwise ldA>=k.
   /// @param B Pointer to real matrix.
   /// @param ldB Column length of the matrix B. If trans = 'N' or 'n' ldB>=n, otherwise ldB>=k.
   /// @param beta Real scalar.
   /// @param C Pointer to real symmetric n-by-n matrix C.
   /// Only the upper or lower triangular part of C is referenced, depending on the value of uplo above.
   /// @param ldC Column length of the matrix C.  ldC>=n
   /// @ingroup BLAS

   template <typename real_t>
   int SYR2K(char uplo, char trans, int_t n, int_t k, real_t alpha, real_t *A, int_t ldA, real_t *B, int_t ldB, real_t beta, real_t *C, int_t ldC)
   {
      using std::toupper;

      const real_t zero(0.0);
      const real_t one(1.0);
      int_t i,j,l;
      real_t *a,*b,*c,*at,*bt;
      real_t s,t;

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
      else if(ldB<((trans=='N')?n:k))
         return -9;
      else if(ldC<n)
         return -12;
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
               b=B;
               for(l=0;l<k;l++)
               {
                  s=alpha*b[j];
                  t=alpha*a[j];
                  for(i=0;i<=j;i++)
                     c[i]+=s*a[i]+t*b[i];
                  a+=ldA;
                  b+=ldB;
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
               b=B;
               for(l=0;l<k;l++)
               {
                  s=alpha*b[j];
                  t=alpha*a[j];
                  for(i=j;i<n;i++)
                     c[i]+=s*a[i]+t*b[i];
                  a+=ldA;
                  b+=ldB;
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
            bt=B;
            for(j=0;j<n;j++)
            {
               a=A;
               b=B;
               for(i=0;i<=j;i++)
               {
                  s=zero;
                  t=zero;
                  for(l=0;l<k;l++)
                  {
                     s+=b[l]*at[l];
                     t+=a[l]*bt[l];
                  }
                  c[i]=alpha*t+alpha*s+beta*c[i];
                  a+=ldA;
                  b+=ldB;
               }
               at+=ldA;
               bt+=ldB;
               c+=ldC;
            }
         }
         else
         {
            c=C;
            at=A;
            bt=B;
            for(j=0;j<n;j++)
            {
               a=A+j*ldA;
               b=B+j*ldB;
               for(i=j;i<n;i++)
               {
                  s=zero;
                  t=zero;
                  for(l=0;l<k;l++)
                  {
                     s+=b[l]*at[l];
                     t+=a[l]*bt[l];
                  }
                  c[i]=alpha*t+alpha*s+beta*c[i];
                  a+=ldA;
                  b+=ldB;
               }
               at+=ldA;
               bt+=ldB;
               c+=ldC;
            }
         }
      }
      return 0;
   }

   /// @brief Performs multiplcation of complex matrices resulting in a symmetric matrix.
   ///
   /// For complex matrices A and B, complex Hermitian matrix C, and complex scalars alpha and beta
   ///
   ///        C := alpha*A*B.'+alpha*B*A.'+beta*C
   /// or
   ///
   ///        C := alpha*A.'*B+alpha*B.'*A+beta*C
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
   ///        if trans = 'N' or 'n' then C := alpha*A*B.'+alpha*B*A.'+beta*C
   ///        if trans = 'T' or 'T' then C := alpha*A.'*B+alpha*B.'*A+beta*C
   /// @param n Specifies the order of the complex symmetric matrix C.  n>=0
   /// @param k Specifies the other dimension of the complex matrices A and B.
   ///
   ///        if trans = 'N' or 'n' then A and B are n-by-k
   ///        if trans = 'T' or 't' then A and B are k-by-n
   /// @param alpha Complex scalar.
   /// @param A Pointer to complex matrix.
   /// @param ldA Column length of the matrix A. If trans = 'N' or 'n' ldA>=n, otherwise ldA>=k.
   /// @param B Pointer to complex matrix.
   /// @param ldB Column length of the matrix B. If trans = 'N' or 'n' ldB>=n, otherwise ldB>=k.
   /// @param beta Complex scalar.
   /// @param C Pointer to complex symmetric n-by-n matrix C.
   /// Only the upper or lower triangular part of C is referenced, depending on the value of uplo above.
   /// @param ldC Column length of the matrix C.  ldC>=n
   /// @ingroup BLAS

   template <typename real_t>
   int SYR2K(char uplo, char trans, int_t n, int_t k, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *B, int_t ldB, complex<real_t> beta, complex<real_t> *C, int_t ldC)
   {
      using std::toupper;

      const complex<real_t> zero(0.0,0.0);
      const complex<real_t> one(1.0,0.0);
      int_t i,j,l;
      complex<real_t> *a,*b,*c,*at,*bt;
      complex<real_t> s,t;

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
      else if(ldB<((trans=='N')?n:k))
         return -9;
      else if(ldC<n)
         return -12;
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
               b=B;
               for(l=0;l<k;l++)
               {
                  s=alpha*b[j];
                  t=alpha*a[j];
                  for(i=0;i<=j;i++)
                     c[i]+=s*a[i]+t*b[i];
                  a+=ldA;
                  b+=ldB;
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
               b=B;
               for(l=0;l<k;l++)
               {
                  s=alpha*b[j];
                  t=alpha*a[j];
                  for(i=j;i<n;i++)
                     c[i]+=s*a[i]+t*b[i];
                  a+=ldA;
                  b+=ldB;
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
            bt=B;
            for(j=0;j<n;j++)
            {
               a=A;
               b=B;
               for(i=0;i<=j;i++)
               {
                  s=zero;
                  t=zero;
                  for(l=0;l<k;l++)
                  {
                     s+=b[l]*at[l];
                     t+=a[l]*bt[l];
                  }
                  c[i]=alpha*t+alpha*s+beta*c[i];
                  a+=ldA;
                  b+=ldB;
               }
               at+=ldA;
               bt+=ldB;
               c+=ldC;
            }
         }
         else
         {
            c=C;
            at=A;
            bt=B;
            for(j=0;j<n;j++)
            {
               a=A+j*ldA;
               b=B+j*ldB;
               for(i=j;i<n;i++)
               {
                  s=zero;
                  t=zero;
                  for(l=0;l<k;l++)
                  {
                     s+=b[l]*at[l];
                     t+=a[l]*bt[l];
                  }
                  c[i]=alpha*t+alpha*s+beta*c[i];
                  a+=ldA;
                  b+=ldB;
               }
               at+=ldA;
               bt+=ldB;
               c+=ldC;
            }
         }
      }
      return 0;
   }
}

#endif

