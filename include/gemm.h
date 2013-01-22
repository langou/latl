//
//  gemm.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/30/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _gemm_h
#define _gemm_h

/// @file gemm.h Computes general matrix-matrix products.

#include <cctype>
#include "latl.h"

namespace latl
{
   /// @brief Computes products of real matrices.
   ///
   /// For real matrices A,B,C and real scalars alpha, beta
   ///
   ///        C := alpha * op(A) * op(B) + beta * C
   ///
   /// is calculated, where op(X) = X or X'.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param transA Specifies whether op(A)=A or A' as follows:
   ///
   ///        if transA = 'N' or 'n' then op(A)=A
   ///        if transA = 'T' or 't' then op(A)=A'
   ///        if transA = 'C' or 'c' then op(A)=A'
   ///
   /// @param transB Specifies whether op(B)=B or B' as follows:
   ///
   ///        if transB = 'N' or 'n' then op(B)=B
   ///        if transB = 'T' or 't' then op(B)=B'
   ///        if transB = 'C' or 'c' then op(B)=B'
   ///
   /// @param m The number of rows of the matrices C and op(A).  m>=0
   /// @param n The number of columns of the matrices C and op(B).  n>=0
   /// @param k The number of columns of the matrix op(A), and the number of rows of the matrix op(B).  k>=0
   /// @param alpha Real scalar.
   /// @param A Pointer to real matrix A, where op(A) is m-by-k.
   /// @param ldA Column length of the matrix A.  If transA='N' or 'n', ldA>=m; otherwise, ldA>=k.
   /// @param B Pointer to real matrix B, where op(B) is k-by-n.
   /// @param ldB Column length of the matrix B.  If transB='N' or 'n', ldB>=k; otherwise, ldA>=n.
   /// @param beta Real scalar.
   /// @param C Pointer to real matrix C, where C is m-by-n.
   /// @param ldC Column length of the matrix C.  ldC>=m.
   /// @ingroup MATM

   template<typename real_t>
   int gemm(char transA, char transB, int_t m, int_t n, int_t k, real_t alpha, real_t *A, int_t ldA, real_t *B, int_t ldB, real_t beta, real_t *C, int_t ldC)
   {
      using std::toupper;
      const real_t zero(0.0);
      const real_t one(1.0);
      int_t i,j,l;
      real_t s,t;
      real_t *a,*b,*c;

      transA=toupper(transA);
      transB=toupper(transB);
      
      if((transA!='N')&&(transA!='T')&&(transA!='C'))
         return -1;
      else if((transB!='N')&&(transB!='T')&&(transB!='C'))
         return -2;
      else if(m<0)
         return -3;
      else if(n<0)
         return -4;
      else if(k<0)
         return -5;
      else if(ldA<((transA!='N')?k:m))
         return -8;
      else if(ldB<((transB!='N')?n:k))
         return -10;
      else if(ldC<m)
         return -13;
      else if((m==0)||(n==0)||(((alpha==zero)||(k==0))&&(beta==one)))
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
      else if((transA=='N')&&(transB=='N'))
      {
         c=C;
         b=B;
         for(j=0;j<n;j++)
         {
            for(i=0;i<m;i++)
               c[i]*=beta;
            a=A;
            for(l=0;l<k;l++)
            {
               for(i=0;i<m;i++)
               {
                  c[i]+=alpha*b[l]*a[i];
               }
               a+=ldA;
            }
            b+=ldB;
            c+=ldC;
         }
      }
      else if((transA!='N')&&(transB=='N'))
      {
         b=B;
         c=C;
         for(j=0;j<n;j++)
         {
            a=A;
            for(i=0;i<m;i++)
            {
               s=zero;
               for(l=0;l<k;l++)
               {
                  s+=a[l]*b[l];
               }
               c[i]=alpha*s+beta*c[i];
               a+=ldA;
            }
            b+=ldB;
            c+=ldC;
         }
      }
      else if((transA=='N')&&(transB!='N'))
      {
         c=C;
         for(j=0;j<n;j++)
         {
            for(i=0;i<m;i++)
               c[i]*=beta;
            a=A;
            b=B;
            for(l=0;l<k;l++)
            {
               t=alpha*b[j];
               for(i=0;i<m;i++)
               {
                  c[i]+=t*a[i];
               }
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
            a=A;
            for(i=0;i<m;i++)
            {
               s=zero;
               b=B;
               for(l=0;l<k;l++)
               {
                  s+=a[l]*b[j];
                  b+=ldB;
               }
               c[i]=alpha*s+beta*c[i];
               a+=ldA;
            }
            c+=ldC;
         }
         
      }
      return 0;
   }
   
   /// @brief Computes products of complex matrices.
   ///
   /// For complex matrices A,B,C and complex scalars alpha, beta
   ///
   ///        C := alpha * op(A) * op(B) + beta * C
   ///
   /// is calculated, where op(X) = X, X' or X.'.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param transA Specifies whether op(A)=A, A.' or A' as follows:
   ///
   ///        if transA = 'N' or 'n' then op(A)=A
   ///        if transA = 'T' or 't' then op(A)=A.'
   ///        if transA = 'C' or 'c' then op(A)=A'
   ///
   /// @param transB Specifies whether op(B)=B, B.' or B' as follows:
   ///
   ///        if transB = 'N' or 'n' then op(B)=B
   ///        if transB = 'T' or 't' then op(B)=B.'
   ///        if transB = 'C' or 'c' then op(B)=B'
   ///
   /// @param m The number of rows of the matrices C and op(A).  m>=0
   /// @param n The number of columns of the matrices C and op(B).  n>=0
   /// @param k The number of columns of the matrix op(A), and the number of rows of the matrix op(B).  k>=0
   /// @param alpha Complex scalar.
   /// @param A Pointer to complex matrix A, where op(A) is m-by-k.
   /// @param ldA Column length of the matrix A.  If transA='N' or 'n', ldA>=m; otherwise, ldA>=k.
   /// @param B Pointer to complex matrix B, where op(B) is k-by-n.
   /// @param ldB Column length of the matrix B.  If transB='N' or 'n', ldB>=k; otherwise, ldA>=n.
   /// @param beta Complex scalar.
   /// @param C Pointer to complex matrix C, where C is m-by-n.
   /// @param ldC Column length of the matrix C.  ldC>=m.
   /// @ingroup MATM

   template<typename real_t>
   int gemm(char transA, char transB, int_t m, int_t n, int_t k, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *B, int_t ldB, complex<real_t> beta, complex<real_t> *C, int_t ldC)
   {
      using std::conj;
      using std::toupper;
      const complex<real_t> zero(0.0,0.0);
      const complex<real_t> one(1.0,0.0);
      int_t i,j,l;
      complex<real_t> *a,*b,*c;
      complex<real_t> s,t;
      
      transA=toupper(transA);
      transB=toupper(transB);

      if((transA!='N')&&(transA!='T')&&(transA!='C'))
         return -1;
      else if((transB!='N')&&(transB!='T')&&(transB!='C'))
         return -2;
      else if(m<0)
         return -3;
      else if(n<0)
         return -4;
      else if(k<0)
         return -5;
      else if(ldA<((transA!='N')?k:m))
         return -8;
      else if(ldB<((transB!='N')?n:k))
         return -10;
      else if(ldC<m)
         return -13;
      else if((m==0)||(n==0)||(((alpha==zero)||(k==0))&&(beta==one)))
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
      else if((transA=='N')&&(transB=='N'))
      {
         c=C;
         b=B;
         for(j=0;j<n;j++)
         {
            for(i=0;i<m;i++)
               c[i]*=beta;
            a=A;
            for(l=0;l<k;l++)
            {
               for(i=0;i<m;i++)
               {
                  c[i]+=alpha*b[l]*a[i];
               }
               a+=ldA;
            }
            b+=ldB;
            c+=ldC;
         }
      }
      else if((transA=='T')&&(transB=='N'))
      {
         b=B;
         c=C;
         for(j=0;j<n;j++)
         {
            a=A;
            for(i=0;i<m;i++)
            {
               s=zero;
               for(l=0;l<k;l++)
               {
                  s+=a[l]*b[l];
               }
               c[i]=alpha*s+beta*c[i];
               a+=ldA;
            }
            b+=ldB;
            c+=ldC;
         }
      }
      else if((transA=='C')&&(transB=='N'))
      {
         b=B;
         c=C;
         for(j=0;j<n;j++)
         {
            a=A;
            for(i=0;i<m;i++)
            {
               s=zero;
               for(l=0;l<k;l++)
               {
                  s+=conj(a[l])*b[l];
               }
               c[i]=alpha*s+beta*c[i];
               a+=ldA;
            }
            b+=ldB;
            c+=ldC;
         }
      }
      else if((transA=='N')&&(transB=='T'))
      {
         c=C;
         for(j=0;j<n;j++)
         {
            for(i=0;i<m;i++)
               c[i]*=beta;
            a=A;
            b=B;
            for(l=0;l<k;l++)
            {
               t=alpha*b[j];
               for(i=0;i<m;i++)
               {
                  c[i]+=t*a[i];
               }
               a+=ldA;
               b+=ldB;
            }
            c+=ldC;
         }
      }
      else if((transA=='N')&&(transB=='C'))
      {
         c=C;
         for(j=0;j<n;j++)
         {
            for(i=0;i<m;i++)
               c[i]*=beta;
            a=A;
            b=B;
            for(l=0;l<k;l++)
            {
               t=alpha*conj(b[j]);
               for(i=0;i<m;i++)
               {
                  c[i]+=t*a[i];
               }
               a+=ldA;
               b+=ldB;
            }
            c+=ldC;
         }
      }
      else if((transA=='T')&&(transB=='T'))
      {
         c=C;
         for(j=0;j<n;j++)
         {
            a=A;
            for(i=0;i<m;i++)
            {
               s=zero;
               b=B;
               for(l=0;l<k;l++)
               {
                  s+=a[l]*b[j];
                  b+=ldB;
               }
               c[i]=alpha*s+beta*c[i];
               a+=ldA;
            }
            c+=ldC;
         }
      }
      else if((transA=='C')&&(transB=='T'))
      {
         c=C;
         for(j=0;j<n;j++)
         {
            a=A;
            for(i=0;i<m;i++)
            {
               s=zero;
               b=B;
               for(l=0;l<k;l++)
               {
                  s+=conj(a[l])*b[j];
                  b+=ldB;
               }
               c[i]=alpha*s+beta*c[i];
               a+=ldA;
            }
            c+=ldC;
         }
      }
      else if((transA=='T')&&(transB=='C'))
      {
         c=C;
         for(j=0;j<n;j++)
         {
            a=A;
            for(i=0;i<m;i++)
            {
               s=zero;
               b=B;
               for(l=0;l<k;l++)
               {
                  s+=a[l]*conj(b[j]);
                  b+=ldB;
               }
               c[i]=alpha*s+beta*c[i];
               a+=ldA;
            }
            c+=ldC;
         }
      }
      else if((transA=='C')&&(transB=='C'))
      {
         c=C;
         for(j=0;j<n;j++)
         {
            a=A;
            for(i=0;i<m;i++)
            {
               s=zero;
               b=B;
               for(l=0;l<k;l++)
               {
                  s+=conj(a[l])*conj(b[j]);
                  b+=ldB;
               }
               c[i]=alpha*s+beta*c[i];
               a+=ldA;
            }
            c+=ldC;
         }
      }
      return 0;
   }
}
#endif
