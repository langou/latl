//
//  geqrt2.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 7/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _geqrt2_h
#define _geqrt2_h

/// @file geqrt2.h Computes a QR factorization of a real matrix, using the compact WY representation of Q.


#include <algorithm>
#include "larfg.h"
#include "gemv.h"
#include "ger.h"
#include "trmv.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes a QR factorization of a real matrix, using the compact WY representation of Q.
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
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param m Specifies the number of rows of the matrix A.  m>=0
   /// @param n Specifies the number of columns of the matrix A.  n>=0
   /// @param A Pointer to real m-by-n matrix A [in,out].  On entry, contains the real m-by-n matrix A.
   /// On exit, the elements on and above the diagonal contain the n-by-n upper triangular matrix R;
   /// the elements below the diagonal are the columns of V, as shown above.
   /// @param ldA Column length of the matrix A.
   /// @param T Pointer to the n-by-n upper triangular factor of the block reflector [out].
   /// The elements on and above the diagonal contain the block reflector T, and the first column of T contains tau.
   /// @param ldT Column length of the matrix T.
   /// @ingroup COMP
   
   template<typename real_t>
   int GEQRT2(int_t m,int_t n,real_t *A,int_t ldA,real_t *T,int_t ldT)
   {
      const real_t zero=0.0;
      const real_t one=1.0;
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;
      else if(ldT<n)
         return -6;
      
      int_t k=std::min(m,n);
      real_t *a=A;
      real_t *b=A+ldA;
      real_t *t=T+(n-1)*ldT;
      for(int_t i=0;i<k;i++)
      {
         LATL::LARFG(m-i,a[i],&(a[i+1]),1,T[i]);
         if(i<n-1)
         {
            real_t temp=a[i];
            a[i]=one;
            LATL::GEMV('T',m-i,n-i-1,one,b+i,ldA,a+i,1,zero,t,1);
            LATL::GER(m-i,n-i-1,-T[i],a+i,1,t,1,b+i,ldA);
            a[i]=temp;
         }
         a+=ldA;
         b+=ldA;
      }
      a=A;
      t=T;
      for(int_t i=1;i<n;i++)
      {
         a+=ldA;
         t+=ldT;
         real_t temp=a[i];
         a[i]=one;
         LATL::GEMV('T',m-i,i,-T[i],A+i,ldA,a+i,1,zero,t,1);
         a[i]=temp;
         LATL::TRMV('U','N','N',i,T,ldT,t,1);
         t[i]=T[i];
         T[i]=zero;
      }
      return 0;
   }

   /// @brief Computes a QR factorization of a complex matrix, using the compact WY representation of Q.
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
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param m Specifies the number of rows of the matrix A.  m>=0
   /// @param n Specifies the number of columns of the matrix A.  n>=0
   /// @param A Pointer to complex m-by-n matrix A [in,out].  On entry, contains the real m-by-n matrix A.
   /// On exit, the elements on and above the diagonal contain the n-by-n upper triangular matrix R;
   /// the elements below the diagonal are the columns of V, as shown above.
   /// @param ldA Column length of the matrix A.
   /// @param T Pointer to the n-by-n upper triangular factor of the block reflector [out].
   /// The elements on and above the diagonal contain the block reflector T, and the first column of T contains tau.
   /// @param ldT Column length of the matrix T.
   /// @ingroup COMP
   
   template<typename real_t>
   int GEQRT2(int_t m,int_t n,complex<real_t> *A,int_t ldA,complex<real_t> *T,int_t ldT)
   {
      using std::conj;
      const complex<real_t> zero(0.0,0.0);
      const complex<real_t> one(1.0,0.0);
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;
      else if(ldT<n)
         return -6;
      
      int_t k=std::min(m,n);
      real_t *a=A;
      real_t *b=A+ldA;
      real_t *t=T+(n-1)*ldT;
      for(int_t i=0;i<k;i++)
      {
         LATL::LARFG(m-i,a[i],&(a[i+1]),1,T[i]);
         if(i<n-1)
         {
            real_t temp=a[i];
            a[i]=one;
            LATL::GEMV('C',m-i,n-i-1,one,b+i,ldA,a+i,ldA,1,zero,t,1);
            LATL::GERC(m-i,n-i-1,-conj(T[i]),a+i,1,t,1,b+i,ldA);
            a[i]=temp;
         }
         a+=ldA;
         b+=ldA;
      }
      a=A;
      t=T;
      for(int_t i=1;i<n;i++)
      {
         a+=ldA;
         t+=ldT;
         real_t temp=a[i];
         a[i]=one;
         LATL::GEMV('C',m-i,i,-T[i],A+i,ldA,a+i,1,zero,t,1);
         a[i]=temp;
         LATL::TRMV('U','N','N',i,T,ldT,t,1);
         t[i]=T[i];
         T[i]=zero;
      }
      return 0;
   }
}
#endif
