//
//  getri.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 1/23/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _getri_h
#define _getri_h

/// @file getri.h Computes the inverse of a matrix using the LU factorization computed by GETRF.

#include "latl.h"
#include "trtri.h"
#include "gemv.h"
#include "gemm.h"
#include "trsm.h"

namespace latl
{
   /// @brief Computes the inverse of a matrix using the LU factorization
   /// computed by GETRF.
   ///
   /// This method inverts U and then computes inv(A) by solving the system
   ///
   ///         inv(A)*L = inv(U) for inv(A)
   ///
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return 1 if the matrix is singular.
   /// @param n Order of the matrix A.  n >= 0.
   /// @param A Real matrix of order n.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param ipiv is integer array with dimension n.  Details of the
   /// interchanges and the block structure of D as determined by GETRF.
   /// @param nb Block size (optional).
   /// @ingroup MATM

   template <typename real_t>
      int_t getri(int_t n, real_t *A, int_t ldA, int_t *ipiv, int_t nb=32)
      {
         const real_t zero(0.0);
         const real_t one(1.0);

         if(n<0)
            return -1;
         else if(ldA<n)
            return -3;
         else if(n==0)
            return 0;

         real_t *work=new real_t[n*nb];

         int_t info=trtri<real_t>('U','N', n, A, ldA);
         if(info != 0)
            return info;

         if((nb<2) || (nb>=n))
         {
            real_t *Aj = A+(n-1)*ldA;
            for(int_t j=n-1;j>=0;j--)
            {
               for(int_t i=j+1;i<n;i++)
               {
                  work[i]=Aj[i];
                  Aj[i]=zero;
               }
               if(j<n-1)
                  gemv<real_t>('N', n, n-j-1, -one, Aj+ldA, ldA, work+j+1, 1, one, Aj, 1);
               Aj -= ldA;
            }
         }
         else
         {
            int_t nn = ( (n-1) / nb ) * nb;
            real_t *Aj = A+nn*ldA;
            for(int_t j=nn+1;j>=1;j-=nb)
            {
               int_t jb = std::min( nb, n-j+1 );
               for(int_t jj=j; jj<=j+jb-1; jj++)
               {
                  real_t *Ajjm1 = A+(jj-1)*ldA;
                  real_t *workjjmj = work+(jj-j)*n;
                  for(int_t i=jj;i<n;i++)
                  {
                     workjjmj[i] = Ajjm1[i];
                     Ajjm1[i] = zero;
                  }
                  Ajjm1 += ldA;
                  workjjmj += n;
               }
               if( j+jb <= n)
                  gemm<real_t>('N', 'N', n, jb, n-j-jb+1, -one, Aj+jb*ldA, ldA, work+j+jb-1, n, one, Aj, ldA);

               trsm<real_t>('R', 'L', 'N', 'U', n, jb, one, work+(j-1), n, Aj, ldA);
               Aj -= nb*ldA;
            }
         }

         real_t *Aj = A+(n-2)*ldA;
         for(int_t j=n-2;j>=0;j--)
         {
            int_t jp = ipiv[j];
            if(jp != j)
               swap<real_t>(n,Aj,1,A+jp*ldA,1);
            Aj -= ldA;
         }

         delete [] work;
         return 0;
      }

   /// @brief Computes the inverse of a matrix using the LU factorization
   /// computed by GETRF.
   ///
   /// This method inverts U and then computes inv(A) by solving the system
   ///
   ///         inv(A)*L = inv(U) for inv(A)
   ///
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return 1 if the matrix is singular.
   /// @param n Order of the matrix A.  n >= 0.
   /// @param A Complex matrix of order n.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param ipiv is integer array with dimension n.  Details of the
   /// interchanges and the block structure of D as determined by GETRF.
   /// @param nb Block size (optional).
   /// @ingroup MATM

   template< typename real_t >
      int_t getri(int_t n, complex<real_t> *A, int_t ldA, int_t *ipiv, int_t nb=32)
      {
         return latl::getri< complex<real_t> > (n, A, ldA, ipiv, nb);
      }

}

#endif

