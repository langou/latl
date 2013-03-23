//
//  herk.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/6/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _herk_h
#define _herk_h

/// @file herk.h Performs complex matrix-matrix multiplication operations resulting in Hermitian matrices.

#include "latl.h"

namespace LATL
{
   /// @brief Performs multiplcation of a complex matrix with its conjugate transpose.
   ///
   /// For complex matrix A, complex Hermitian matrix C, and real scalars alpha and beta
   ///
   ///        C := alpha*A*A'+beta*C  or  C := alpha*A'*A+beta*C
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the Hermitian matrix C
   /// is to be referenced:  
   ///
   ///        if uplo = 'U' or 'u' then C is upper triangular,
   ///        if uplo = 'L' or 'l' then C is lower triangular.
   /// @param trans Specifies the operation to be perfomed as follows:
   ///
   ///        if trans = 'N' or 'n' then C := alpha*A*A'+beta*C
   ///        if trans = 'C' or 'c' then C := alpha*A'*A+beta*C
   /// @param n Specifies the order of the complex Hermitian matrix C.  n>=0
   /// @param k Specifies the other dimension of the complex matrix A (see below). k>=0
   /// @param alpha Real scalar.
   /// @param A Pointer to complex matrix.
   ///
   ///        if trans = 'N' or 'n' then A is n-by-k
   ///        if trans = 'C' or 'c' then A is k-by-n
   /// @param ldA Column length of the matrix A. If trans = 'N' or 'n' ldA>=n, otherwise ldA>=k.
   /// @param beta Real scalar.
   /// @param C Pointer to complex Hermitian n-by-n matrix C.  
   /// Only the upper or lower triangular part of C is referenced, depending on the value of uplo above.
   /// @param ldC Column length of the matrix C.  ldC>=n
   /// @ingroup MATM

   template <typename real_t>
   int HERK(char uplo, char trans, int_t n, int_t k, real_t alpha, complex<real_t> *A, int_t ldA, real_t beta, complex<real_t> *C, int_t ldC)
   {
      using std::real;
      using std::conj;
      const complex<real_t> zero(0.0,0.0);
      const real_t rzero=0.0;
      const real_t one=1.0;
      int_t i,j,l;
      complex<real_t> *a,*c,*at;
      complex<real_t> t;
      real_t s;
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='C'))
         return -2;
      else if(n<0)
         return -3;
      else if(k<0)
         return -4;
      else if(ldA<((trans=='N')?n:k))
         return -7;
      else if(ldC<n)
         return -10;
      else if((n==0)||(((alpha==rzero)||(k==0))&&(beta==one)))
         return 0;
      
      if(alpha==zero)
      {
         if(uplo=='U')
         {
            c=C;
            for(j=0;j<n;j++)
            {
               for(i=0;i<j;i++)
                  c[i]*=beta;
               c[j]=beta*real(c[j]);
               c+=ldC;
            }
            
         }
         else
         {
            c=C;
            for(j=0;j<n;j++)
            {
               c[j]=beta*real(c[j]);
               for(i=j+1;i<n;i++)
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
               for(i=0;i<j;i++)
                  c[i]*=beta;
               c[j]=beta*real(c[j]);
               a=A;
               for(l=0;l<k;l++)
               {
                  t=alpha*conj(a[j]);
                  for(i=0;i<j;i++)
                     c[i]+=t*a[i];
                  c[j]=real(c[j])+real(t*a[i]);
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
               c[j]=beta*real(c[j]);
               for(i=j+1;i<n;i++)
                  c[i]*=beta;
               a=A;
               for(l=0;l<k;l++)
               {
                  t=alpha*conj(a[j]);
                  c[j]=real(c[j])+real(t*a[j]);
                  for(i=j+1;i<n;i++)
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
               for(i=0;i<j;i++)
               {
                  t=zero;
                  for(l=0;l<k;l++)
                     t+=conj(a[l])*at[l];
                  c[i]=alpha*t+beta*c[i];
                  a+=ldA;
               }
               s=rzero;
               for(l=0;l<k;l++)
                  s+=real(conj(at[l])*at[l]);
               c[j]=alpha*s+beta*real(c[j]);
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
               s=rzero;
               for(l=0;l<k;l++)
                  s+=real(conj(at[l])*at[l]);
               c[j]=alpha*s+beta*real(c[j]);
               a=A+(j+1)*ldA;
               for(i=j+1;i<n;i++)
               {
                  t=zero;
                  for(l=0;l<k;l++)
                     t+=conj(a[l])*at[l];
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
   
   /// @brief Performs multiplcation of complex matrices resulting in a Hermitian matrix.
   ///
   /// For complex matrices A and B, complex Hermitian matrix C, and complex scalar alpha and real 
   /// scalar beta
   ///
   ///        C := alpha*A*B'+conj(alpha)*B*A'+beta*C
   /// or
   ///
   ///        C := alpha*A'*B+conj(alpha)*B'*A+beta*C
   /// is computed.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the Hermitian matrix C
   /// is to be referenced:  
   ///
   ///        if uplo = 'U' or 'u' then C is upper triangular,
   ///        if uplo = 'L' or 'l' then C is lower triangular.
   /// @param trans Specifies the operation to be perfomed as follows:
   ///
   ///        if trans = 'N' or 'n' then C := alpha*A*B'+conj(alpha)*B*A'+beta*C
   ///        if trans = 'C' or 'c' then C := alpha*A'*B+conj(alpha)*B'*A+beta*C
   /// @param n Specifies the order of the complex Hermitian matrix C.  n>=0
   /// @param k Specifies the other dimension of the complex matrices A and B.
   ///
   ///        if trans = 'N' or 'n' then A and B are n-by-k
   ///        if trans = 'C' or 'c' then A and B are k-by-n
   /// @param alpha Complex scalar.
   /// @param A Pointer to complex matrix.
   /// @param ldA Column length of the matrix A. If trans = 'N' or 'n' ldA>=n, otherwise ldA>=k.
   /// @param B Pointer to complex matrix.
   /// @param ldB Column length of the matrix B. If trans = 'N' or 'n' ldB>=n, otherwise ldB>=k.
   /// @param beta Real scalar.
   /// @param C Pointer to complex Hermitian n-by-n matrix C.  
   /// Only the upper or lower triangular part of C is referenced, depending on the value of uplo above.
   /// @param ldC Column length of the matrix C.  ldC>=n
   /// @ingroup MATM

   template <typename real_t>
   int HERK(char uplo, char trans, int_t n, int_t k, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *B, int_t ldB, complex<real_t> beta, complex<real_t> *C, int_t ldC)
   {
      using std::real;
      using std::conj;
      const complex<real_t> zero(0.0,0.0);
      const complex<real_t> one(1.0,0.0);
      int_t i,j,l;
      complex<real_t> *a,*b,*c,*at,*bt;
      complex<real_t> s,t;
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='C'))
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
               for(i=0;i<j;i++)
                  c[i]*=beta;
               c[j]=beta*real(c[j]);
               c+=ldC;
            }
            
         }
         else
         {
            c=C;
            for(j=0;j<n;j++)
            {
               c[j]=beta*real(c[j]);
               for(i=j+1;i<n;i++)
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
               for(i=0;i<j;i++)
                  c[i]*=beta;
               c[j]=beta*real(c[j]);
               a=A;
               b=B;
               for(l=0;l<k;l++)
               {
                  s=alpha*conj(b[j]);
                  t=conj(alpha*a[j]);
                  for(i=0;i<j;i++)
                     c[i]+=s*a[i]+t*b[i];
                  c[j]=real(c[j])+real(a[j]*s+b[j]*t);
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
               c[j]=beta*real(c[j]);
               for(i=j+1;i<n;i++)
                  c[i]*=beta;
               a=A;
               b=B;
               for(l=0;l<k;l++)
               {
                  s=alpha*conj(b[j]);
                  t=conj(alpha*a[j]);
                  c[j]=real(c[j])+real(s*a[j]+t*b[j]);
                  for(i=j+1;i<n;i++)
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
                     s+=conj(b[l])*at[l];
                     t+=conj(a[l])*bt[l];
                  }
                  if(i<j)
                     c[i]=alpha*t+conj(alpha)*s+beta*c[i];
                  else
                     c[j]=real(alpha*t+conj(alpha)*s)+beta*real(c[j]);
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
                     s+=conj(b[l])*at[l];
                     t+=conj(a[l])*bt[l];                  
                  }
                  if(i>j)
                     c[i]=alpha*t+conj(alpha)*s+beta*c[i];
                  else
                     c[j]=real(alpha*t+conj(alpha)*s)+beta*real(c[j]);
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
