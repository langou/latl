//
//  geqrt.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 2/25/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _geqrt_h
#define _geqrt_h

/// @file geqrt.h Computes a QR factoriztion of a rectangular matrix.

#include <algorithm>
#include "geqrt3.h"
#include "larfb.h"
#include "latl.h"

namespace LATL
{

   /// @brief Computes a QR factorization of a real matrix using blocked algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param m Number of rows of matrix A.  m>=0
   /// @param n Number of columns of matrix A. n>=0
   /// @param[in,out] A Real m-by-n matrix.  On exit, the portion above the diagonal contains the triangular
   /// factor R, and the portion below the diagonal contains the Householder vectors stored as column vectors.
   /// @param ldA Leading dimension of the matrix A.  ldA>=m.
   /// @param T Upper triangular part of the block reflectors, stored as nb-by-n matrix.
   /// @param ldT Leading dimension of the matrix T.  ldT>=nb
   /// @param nb Block size.  0<nb<=min(m,n)
   /// @param W Workspace vector of length n*min(m,n) (optional).  If not provided, workspace is managed internally.

   template<typename real_t>
   int GEQRT(int_t m, int_t n, real_t *A, int_t ldA, real_t *T, int_t ldT, int_t nb, real_t *W=NULL)
   {
      using std::min;
      int_t k=min(m,n);
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;
      else if(ldT<nb)
         return -6;
      else if((nb<1)||(nb>k))
         return -7;
      else if(k==0)
         return 0;

      bool allocate=(W==NULL)?1:0;
      if(allocate)
         W=new real_t[k*n];
      real_t *B=A;
      for(int_t i=0;i<k;i+=nb)
      {
         int_t ib=min(k-i,nb);
         B+=ib*ldA;
         GEQRT3(m-i,ib,A+i,ldA,T,ldT);
         if(i+ib<n)
            LARFB('L','T','F','C',m-i,n-i-ib,ib,A+i,ldA,T,ldT,B+i,ldA,W);
         A=B;
      }
      if(allocate)
         delete [] W;
      return 0;
   }

   /// @brief Computes a QR factorization of a complex matrix using blocked algorithm.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @param m Number of rows of matrix A.  m>=0
   /// @param n Number of columns of matrix A. n>=0
   /// @param[in,out] A Complex m-by-n matrix.  On exit, the portion above the diagonal contains the triangular
   /// factor R, and the portion below the diagonal contains the Householder vectors stored as column vectors.
   /// @param ldA Leading dimension of the matrix A.  ldA>=m.
   /// @param T Upper triangular part of the block reflectors, stored as nb-by-n matrix.
   /// @param ldT Leading dimension of the matrix T.  ldT>=nb
   /// @param nb Block size.  0<nb<=min(m,n)
   /// @param W Workspace vector of length n*min(m,n) (optional).  If not provided, workspace is managed internally.

   template<typename real_t>
   int GEQRT(int_t m, int_t n, complex<real_t> *A, int_t ldA, complex<real_t> *T, int_t ldT, int_t nb, complex<real_t> *W=NULL)
   {
      using std::min;
      int_t k=min(m,n);
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;
      else if(ldT<nb)
         return -6;
      else if((nb<1)||(nb>k))
         return -7;
      else if(k==0)
         return 0;

      bool allocate=(W==NULL)?1:0;
      if(allocate)
         W=new complex<real_t>[k*n];
      complex<real_t> *B=A;
      for(int_t i=0;i<k;i+=nb)
      {
         int_t ib=min(k-i,nb);
         B+=ib*ldA;
         GEQRT3(m-i,ib,A+i,ldA,T,ldT);
         if(i+ib<n)
            LARFB('L','C','F','C',m-i,n-i-ib,ib,A+i,ldA,T,ldT,B+i,ldA,W);
         A=B;
      }
      if(allocate)
         delete [] W;
      return 0;
   }
   
}



#endif
