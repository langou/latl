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
   /// @brief Recursively computes a QR factorization of a real matrix, using the compact WY representation of Q. 
   ///
   /// The matrix V stores the elementary reflectors H(i) in the i-th column
   /// below the diagonal. For example, if M=5 and N=3, the matrix V is
   ///
   ///               V = (  1       )
   ///                   ( v1  1    )
   ///                   ( v1 v2  1 )
   ///                   ( v1 v2 v3 )
   ///                   ( v1 v2 v3 )
   /// where the vi's represent the vectors which define H(i), which are returned
   /// in the matrix A.  The 1's along the diagonal of V are not stored in A.  The
   /// block reflector H is then given by
   ///
   ///               H = I - V T V'
   /// For details of the algorithm, see Elmroth and Gustavson, IBM J. Res. Develop. Vol 44 No. 4 July 2000.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param m Specifies the number of rows of the matrix A.  m>=n
   /// @param n Specifies the number of columns of the matrix A.  n>=0
   /// @param A Pointer to real m-by-n matrix A [in,out].  On entry, contains the real m-by-n matrix A.
   /// On exit, the elements on and above the diagonal contain the n-by-n upper triangular matrix R;
   /// the elements below the diagonal are the columns of V, as shown above.
   /// @param ldA Column length of the matrix A.
   /// @param T Pointer to the n-by-n upper triangular factor of the block reflector [out].  
   /// The elements on and above the diagonal contain the block reflector T; the elements below the diagonal are not used.
   /// @param ldT Column length of the matrix T.
   /// @ingroup QRF
   
   template<typename real_t>
   int GEQRT3(int_t m,int_t n,real_t *A,int_t ldA,real_t *T,int_t ldT)
   {
      using LATL::LARFG;
      using LATL::GEMM;
      using LATL::TRMM;
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
   
   /// @brief Recursively computes a QR factorization of a complex matrix, using the compact WY representation of Q. 
   ///
   /// The matrix V stores the elementary reflectors H(i) in the i-th column
   /// below the diagonal. For example, if M=5 and N=3, the matrix V is
   ///
   ///               V = (  1       )
   ///                   ( v1  1    )
   ///                   ( v1 v2  1 )
   ///                   ( v1 v2 v3 )
   ///                   ( v1 v2 v3 )
   /// where the vi's represent the vectors which define H(i), which are returned
   /// in the matrix A.  The 1's along the diagonal of V are not stored in A.  The
   /// block reflector H is then given by
   ///
   ///               H = I - V T V'
   /// For details of the algorithm, see Elmroth and Gustavson, IBM J. Res. Develop. Vol 44 No. 4 July 2000.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param m Specifies the number of rows of the matrix A.  m>=n
   /// @param n Specifies the number of columns of the matrix A.  n>=0
   /// @param A Pointer to complex m-by-n matrix A [in,out].  On entry, contains the complex m-by-n matrix A.
   /// On exit, the elements on and above the diagonal contain the n-by-n upper triangular matrix R;
   /// the elements below the diagonal are the columns of V, as shown above.
   /// @param ldA Column length of the matrix A.
   /// @param T Pointer to the n-by-n upper triangular factor of the block reflector [out].  
   /// The elements on and above the diagonal contain the block reflector T; the elements below the diagonal are not used.
   /// @param ldT Column length of the matrix T.
   /// @ingroup QRF
   
   template<typename real_t>
   int GEQRT3(int_t m,int_t n,complex<real_t> *A,int_t ldA,complex<real_t> *T,int_t ldT)
   {
      using std::conj;
      using LATL::LARFG;
      using LATL::GEMM;
      using LATL::TRMM;
      const complex<real_t> one(1.0,0.0);
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
         complex<real_t> *A2;
         complex<real_t> *T2;
         
         A2=A+n1*ldA;
         T2=T+n1*ldT;
         for(int_t j=0;j<n2;j++)
         {
            for(int_t i=0;i<n1;i++)
               T2[i]=A2[i];
            A2+=ldA;
            T2+=ldT;
         }
         
         TRMM('L','L','C','U',n1,n2,one,A,ldA,T1,ldT);
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
