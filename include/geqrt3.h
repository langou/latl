//
//  geqrt3.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 7/18/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _geqrt3_h
#define _geqrt3_h

/// @file geqrt3.h Recursive QR factorization using the compact WY representation of Q.


#include "larfg.h"
#include "trmm.h"
#include "gemm.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes a QR factorization of a real matrix using recursive algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param m Number of rows of matrix A.  m>=n
   /// @param n Number of columns of matrix A. n>=0
   /// @param[in,out] A Real m-by-n matrix.  On exit, the portion above the diagonal contains the triangular
   /// factor R, and the portion below the diagonal contains the Householder vectors stored as column vectors.
   /// @param ldA Leading dimension of the matrix A.  ldA>=m.
   /// @param T Upper triangular part of the block reflector of order n.
   /// @param ldT Leading dimension of the matrix T.  ldT>=n

   template<typename real_t>
   int GEQRT3(int_t m,int_t n,real_t *A,int_t ldA,real_t *T,int_t ldT)
   {
      const real_t one(1.0);
      if(n<0)
         return -2;
      else if(m<n)
         return -1;
      else if(ldA<m)
         return -4;
      else if(ldT<n)
         return -6;

      if(n==1)
      {
         LARFG(m,A[0],A+1,1,T[0]);
      }
      else
      {
         int_t n1=n/2;
         int_t n2=n-n1;

         GEQRT3(m,n1,A,ldA,T,ldT);

         real_t *A1=A+n1*ldA;
         real_t *T1=T+n1*ldT;
         real_t *A2=A1;
         real_t *T2=T1;

         for(int_t j=0;j<n2;j++)
         {
            for(int_t i=0;i<n1;i++)
               T2[i]=A2[i];
            A2+=ldA;
            T2+=ldT;
         }

         TRMM('L','L','T','U',n1,n2,one,A,ldA,T1,ldT);
         GEMM('T','N',n1,n2,m-n1,one,A+n1,ldA,A1+n1,ldA,one,T1,ldT);
         TRMM('L','U','T','N',n1,n2,one,T,ldT,T1,ldT);
         GEMM('N','N',m-n1,n2,n1,-one,A+n1,ldA,T1,ldT,one,A1+n1,ldA);
         TRMM('L','L','N','U',n1,n2,one,A,ldA,T1,ldT);

         A2=A1;
         T2=T1;
         for(int_t j=0;j<n2;j++)
         {
            for(int_t i=0;i<n1;i++)
               A2[i]-=T2[i];
            A2+=ldA;
            T2+=ldT;
         }

         GEQRT3(m-n1,n2,A1+n1,ldA,T1+n1,ldT);

         A2=A;
         for(int_t i=0;i<n1;i++)
         {
            T2=T1;
            for(int_t j=0;j<n2;j++)
            {
               T2[i]=A2[j+n1];
               T2+=ldT;
            }
            A2+=ldA;
         }

         TRMM('R','L','N','U',n1,n2,one,A1+n1,ldA,T1,ldT);
         GEMM('T','N',n1,n2,m-n,one,A+n,ldA,A1+n,ldA,one,T1,ldT);
         TRMM('L','U','N','N',n1,n2,-one,T,ldT,T1,ldT);
         TRMM('R','U','N','N',n1,n2,one,T1+n1,ldT,T1,ldT);
      }
      return 0;
   }

   /// @brief Computes a QR factorization of a complex matrix using recursive algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param m Number of rows of matrix A.  m>=n
   /// @param n Number of columns of matrix A. n>=0
   /// @param[in,out] A Complex m-by-n matrix.  On exit, the portion above the diagonal contains the triangular
   /// factor R, and the portion below the diagonal contains the Householder vectors stored as column vectors.
   /// @param ldA Leading dimension of the matrix A.  ldA>=m.
   /// @param T Upper triangular part of the block reflector of order n.
   /// @param ldT Leading dimension of the matrix T.  ldT>=n

   template<typename real_t>
   int GEQRT3(int_t m,int_t n,complex<real_t> *A,int_t ldA,complex<real_t> *T,int_t ldT)
   {
      using std::conj;
      const complex<real_t> one(1.0);
      if(n<0)
         return -2;
      else if(m<n)
         return -1;
      else if(ldA<m)
         return -4;
      else if(ldT<n)
         return -6;

      if(n==1)
      {
         LARFG(m,A[0],A+1,1,T[0]);
      }
      else
      {
         int_t n1=n/2;
         int_t n2=n-n1;

         GEQRT3(m,n1,A,ldA,T,ldT);

         complex<real_t> *A1=A+n1*ldA;
         complex<real_t> *T1=T+n1*ldT;
         complex<real_t> *A2=A1;
         complex<real_t> *T2=T1;

         for(int_t j=0;j<n2;j++)
         {
            for(int_t i=0;i<n1;i++)
               T2[i]=A2[i];
            A2+=ldA;
            T2+=ldT;
         }

         TRMM('L','L','T','U',n1,n2,one,A,ldA,T1,ldT);
         GEMM('C','N',n1,n2,m-n1,one,A+n1,ldA,A1+n1,ldA,one,T1,ldT);
         TRMM('L','U','C','N',n1,n2,one,T,ldT,T1,ldT);
         GEMM('N','N',m-n1,n2,n1,-one,A+n1,ldA,T1,ldT,one,A1+n1,ldA);
         TRMM('L','L','N','U',n1,n2,one,A,ldA,T1,ldT);

         A2=A1;
         T2=T1;
         for(int_t j=0;j<n2;j++)
         {
            for(int_t i=0;i<n1;i++)
               A2[i]-=T2[i];
            A2+=ldA;
            T2+=ldT;
         }

         GEQRT3(m-n1,n2,A1+n1,ldA,T1+n1,ldT);

         A2=A;
         for(int_t i=0;i<n1;i++)
         {
            T2=T1;
            for(int_t j=0;j<n2;j++)
            {
               T2[i]=conj(A2[j+n1]);
               T2+=ldT;
            }
            A2+=ldA;
         }

         TRMM('R','L','N','U',n1,n2,one,A1+n1,ldA,T1,ldT);
         GEMM('C','N',n1,n2,m-n,one,A+n,ldA,A1+n,ldA,one,T1,ldT);
         TRMM('L','U','N','N',n1,n2,-one,T,ldT,T1,ldT);
         TRMM('R','U','N','N',n1,n2,one,T1+n1,ldT,T1,ldT);
      }
      return 0;
   }
}

#endif
