//
//  laset.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 7/9/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _laset_h
#define _laset_h

/// @file laset.h Initializes a matrix to diagonal and off-diagonal values.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Initializes a real matrix with diagonal values of beta and off-diagonal values of alpha.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the matrix A is upper or lower triangular:
   ///
   ///        'U': A is assumed to be upper triangular; elements below the diagonal are not referenced.
   ///        'L': A is assumed to be lower triangular; elements above the diagonal are not referenced.
   ///        otherwise, A is assumed to be a full matrix.
   /// @param m The number of rows of the matrix A.
   /// @param n The number of columns of the matrix A.
   /// @param alpha Value to assign to the off-diagonal elements of A.
   /// @param beta Value to assign to the diagonal elements of A.
   /// @param A Pointer to real matrix A [out].
   /// @param ldA Column length of the matrix A.
   /// @ingroup AUX
   
   template<typename real_t>
   int LASET(char uplo,int_t m,int_t n,real_t alpha,real_t beta,real_t *A,int_t ldA)
   {
      using std::toupper;
      uplo=toupper(uplo);
      if(m<0)
         return -2;
      else if(n<0)
         return -3;
      else if(ldA<m)
         return -7;
      
      if(uplo=='U')
      {
         for(int_t j=0;j<n;j++)
         {
            for(int_t i=0;(i<m)&&(i<=j);i++)
               A[i]=(i==j)?beta:alpha;
            A+=ldA;
         }
      }
      else if(uplo=='L')
      {
         for(int_t j=0;(j<n)&&(j<m);j++)
         {
            for(int_t i=j;i<m;i++)
               A[i]=(i==j)?beta:alpha;
            A+=ldA;
         }
      }
      else
      {
         for(int_t j=0;j<n;j++)
         {
            for(int_t i=0;i<m;i++)
               A[i]=(i==j)?beta:alpha;
            A+=ldA;
         }
      }
      return 0;
   }
   
   /// @brief Initializes a complex matrix with diagonal values of beta and off-diagonal values of alpha.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the matrix A is upper or lower triangular:
   ///
   ///        'U': A is assumed to be upper triangular; elements below the diagonal are not referenced.
   ///        'L': A is assumed to be lower triangular; elements above the diagonal are not referenced.
   ///        otherwise, A is assumed to be a full matrix.
   /// @param m The number of rows of the matrix A.
   /// @param n The number of columns of the matrix A.
   /// @param alpha Value to assign to the off-diagonal elements of A.
   /// @param beta Value to assign to the diagonal elements of A.
   /// @param A Pointer to complex matrix A [out].
   /// @param ldA Column length of the matrix A.
   /// @ingroup AUX
   
   template<typename real_t>
   int LASET(char uplo,int_t m,int_t n,complex<real_t> alpha,complex<real_t> beta,complex<real_t> *A,int_t ldA)
   {
      return LASET< complex<real_t> >(uplo,m,n,alpha,beta,A,ldA);
   }
}

#endif
