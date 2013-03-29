//
//  lapmt.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lapmt_h
#define _lapmt_h

/// @file lapmt.h Permutes the columns of a matrix.

#include "latl.h"

namespace LATL
{
   /// @brief Permutes the columns of a real matrix.
   ///
   /// Rearranges the columns of the m-by-n matrix X as specified
   /// by the permutation k[0],k[1],...,k[m-1].
   /// If forw is true, a forward permutation is performed:
   ///
   ///        col(k[i]) -> col(i) for each i=0,...,n-1
   /// If forw is false, a backward permutation is performed:
   ///
   ///        col(i) -> col(k[i]) for each i=0,...,n-1
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param forw If true, a forward permutation is performed; if false a backward permutation is performed.
   /// @param m The number of rows of the matrix X.
   /// @param n The number of columns of the matrix X.
   /// @param X Pointer to the real matrix X. On exit, the columns of X are permuted as specified by k.
   /// @param ldX Column length of the matrix X.  ldX>=m
   /// @param k Pointer to permutation vector of length m.
   /// @ingroup AUX

   template<typename real_t>
   int LAPMT(bool forw,int_t m,int_t n,real_t *X,int_t ldX,int_t *k)
   {
      if(m<0)
         return -2;
      else if(n<0)
         return -3;
      else if(ldX<m)
         return -5;
      
      bool *move=new bool[n];
      for(int_t i=0;i<n;i++)
         move[i]=1;
      
      if(forw)
      {
         for(int_t i=0;i<n;i++)
         {
            if(move[i])
            {
               int_t j=i;
               move[j]=!move[j];
               int_t input=k[j];
               while(move[input])
               {
                  real_t *a=X+j*ldX;
                  real_t *b=X+input*ldX;
                  for(int_t l=0;l<m;l++)
                  {
                     real_t temp=a[l];
                     a[l]=b[l];
                     b[l]=temp;
                  }
                  move[input]=!move[input];
                  j=input;
                  input=k[input];
               }
            }
         }
      }
      else
      {
         for(int_t i=0;i<n;i++)
         {
            if(move[i])
            {
               move[i]=!move[i];
               int_t j=k[i];
               while(j!=i)
               {
                  real_t *a=X+i*ldX;
                  real_t *b=X+j*ldX;
                  for(int_t l=0;l<m;l++)
                  {
                     real_t temp=a[l];
                     a[l]=b[l];
                     b[l]=temp;
                  }
                  move[j]=!move[j];
                  j=k[j];
               }
            }
         }
         
      }
      delete [] move;
      return 0;
   }
   
   /// @brief Permutes the columns of a complex matrix.
   ///
   /// Rearranges the columns of the m-by-n matrix X as specified
   /// by the permutation k[0],k[1],...,k[m-1].
   /// If forw is true, a forward permutation is performed:
   ///
   ///        col(k[i]) -> col(i) for each i=0,...,n-1
   /// If forw is false, a backward permutation is performed:
   ///
   ///        col(i) -> col(k[i]) for each i=0,...,n-1
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param forw If true, a forward permutation is performed; if false a backward permutation is performed.
   /// @param m The number of rows of the matrix X.
   /// @param n The number of columns of the matrix X.
   /// @param X Pointer to the complex matrix X. On exit, the columns of X are permuted as specified by k.
   /// @param ldX Column length of the matrix X.  ldX>=m
   /// @param k Pointer to permutation vector of length m.
   /// @ingroup AUX
   
   template<typename real_t>
   int LAPMT(bool forw,int_t m,int_t n,complex<real_t> *X,int_t ldX,int_t *k)
   {
      return LAPMT< complex<real_t> >(forw,m,n,X,ldX,k);
   }
}

#endif
