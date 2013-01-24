//
//  ttmm.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/23/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _ttmm_h
#define _ttmm_h

/// @file ttmm.h Computes product of triangular matrices.

#include <cctype>
#include "latl.h"

namespace latl
{
   /// @brief Computes products of real triangular matrices.
   ///
   /// For real matrices A,B,C and real scalars alpha, beta
   ///
   ///        C <- alpha * A * B + beta * C
   /// is calculated, where A,B, and C are all either upper or lower triangular.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A,B,C are upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A,B,C are upper triangular
   ///        if uplo = 'L' or 'l' then A,B,C are lower triangular
   /// If uplo is set to 'U' or 'u', the portion of each matrix below the diagonal is not accessed;
   /// if uplo is set to 'L' or 'l', the portion of each matrix above the diagonal is not accessed.
   /// @param n The number order of matrices A,B,C.
   /// @param alpha Real scalar.
   /// @param A Real triangular matrix of order n.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param B Real triangular matrix of order n.
   /// @param ldB Column length of the matrix B.  ldA>=n.
   /// @param beta Real scalar.
   /// @param C Real triangular matrix of order n.  On exit, C contains the result of alpha*A*B+beta*C.
   /// @param ldC Column length of the matrix C.  ldC>=n.
   /// @ingroup MATM
   
   template <typename real_t>
   int ttmm(char uplo, int_t n, real_t alpha, real_t *A, int_t ldA, real_t *B, int_t ldB, real_t beta, real_t *C, int_t ldC)
   {
      const real_t zero=0.0;
      using std::toupper;
      uplo=toupper(uplo);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<n)
         return -5;
      else if(ldB<n)
         return -7;
      else if(ldC<n)
         return -10;
      else if(n==0)
         return 0;
      
      if(uplo=='U')
      {
         if(alpha==zero)
         {
            real_t *c=C;
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=0;i<=j;i++)
                  c[i]*=beta;
               c+=ldC;
            }
         }
         else
         {
            real_t *c=C;
            real_t *b=B;
            for(int_t j=0;j<n;j++)
            {
               real_t *ai=A;
               for(int_t i=0;i<=j;i++)
               {
                  c[i]*=beta;
                  real_t *a=ai;
                  for(int_t k=i;k<=j;k++)
                  {
                     c[i]+=alpha*a[i]*b[k];
                     a+=ldA;
                  }
                  ai+=ldA;
               }
               b+=ldB;
               c+=ldC;
            }
         }
      }
      else
      {
         if(alpha==zero)
         {
            real_t *c=C;
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=j;i<n;i++)
                  c[i]*=beta;
               c+=ldC;
            }
         }
         else
         {
            real_t *c=C;
            real_t *b=B;
            real_t *aj=A;
            for(int_t j=0;j<n;j++)
            {
               for(int_t i=j;i<n;i++)
               {
                  c[i]*=beta;
                  real_t *a=aj;
                  for(int_t k=j;k<=i;k++)
                  {
                     c[i]+=alpha*a[i]*b[k];
                     a+=ldA;
                  }
               }
               aj+=ldA;
               b+=ldB;
               c+=ldC;
            }
         }
      }
      return 0;
   }
   
   /// @brief Computes products of complex triangular matrices.
   ///
   /// For complex matrices A,B,C and real scalars alpha, beta
   ///
   ///        C <- alpha * A * B + beta * C
   /// is calculated, where A,B, and C are all either upper or lower triangular.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A,B,C are upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A,B,C are upper triangular
   ///        if uplo = 'L' or 'l' then A,B,C are lower triangular
   /// If uplo is set to 'U' or 'u', the portion of each matrix below the diagonal is not accessed;
   /// if uplo is set to 'L' or 'l', the portion of each matrix above the diagonal is not accessed.
   /// @param n The number order of matrices A,B,C.
   /// @param alpha Complex scalar.
   /// @param A Complex triangular matrix of order n.
   /// @param ldA Column length of the matrix A.  ldA>=n.
   /// @param B Complex triangular matrix of order n.
   /// @param ldB Column length of the matrix B.  ldA>=n.
   /// @param beta Complex scalar.
   /// @param C Complex triangular matrix of order n.  On exit, C contains the result of alpha*A*B+beta*C.
   /// @param ldC Column length of the matrix C.  ldC>=n.
   /// @ingroup MATM
   
   template <typename real_t>
   int ttmm(char uplo, int_t n, complex<real_t> alpha, complex<real_t> *A, int_t ldA, complex<real_t> *B, int_t ldB, complex<real_t> beta, complex<real_t> *C, int_t ldC)
   {
      return ttmm< complex<real_t> >(uplo,n,alpha,A,ldA,B,ldB,beta,C,ldC);
   }
}

#endif
