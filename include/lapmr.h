//
//  lapmr.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lapmr_h
#define _lapmr_h

/// @file lapmr.h Permutes the rows of a matrix.

#include "latl.h"

namespace LATL
{
   /// @brief Permutes the rows of a real matrix.
   ///
   /// Rearranges the rows of the m-by-n matrix X as specified
   /// by the permutation k[0],k[1],...,k[m-1].
   /// If forw is true, a forward permutation is performed:
   ///
   ///        row(k[i]) -> row(i) for each i=0,...,m-1
   /// If forw is false, a backward permutation is performed:
   ///
   ///        row(i) -> row(k[i]) for each i=0,...,m-1
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param forw If true, a forward permutation is performed; if false a backward permutation is performed.
   /// @param m The number of rows of the matrix X.
   /// @param n The number of columns of the matrix X.
   /// @param X Pointer to the real matrix X. On exit, the rows of X are permuted as specified by k.
   /// @param ldX Column length of the matrix X.  ldX>=m
   /// @param k Pointer to permutation vector of length m.
   /// @ingroup MAT

   template<typename real_t>
   int LAPMR(bool forw,int_t m,int_t n,real_t *X,int_t ldX,int_t *k)
   {
      if(m<0)
         return -2;
      else if(n<0)
         return -3;
      else if(ldX<m)
         return -5;
      
      bool *move=new bool[m];
      for(int_t i=0;i<m;i++)
         move[i]=1;
      
      if(forw)
      {
         for(int_t i=0;i<m;i++)
         {
            if(move[i])
            {
               int_t j=i;
               move[j]=!move[j];
               int_t input=k[j];
               while(move[input])
               {
                  real_t *x=X;
                  for(int_t l=0;l<n;l++)
                  {
                     real_t temp=x[j];
                     x[j]=x[input];
                     x[input]=temp;
                     x+=ldX;
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
         for(int_t i=0;i<m;i++)
         {
            if(move[i])
            {
               move[i]=!move[i];
               int_t j=k[i];
               while(j!=i)
               {
                  real_t *x=X;
                  for(int_t l=0;l<n;l++)
                  {
                     real_t temp=x[i];
                     x[i]=x[j];
                     x[j]=temp;
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
   
   /// @brief Permutes the rows of a complex matrix.
   ///
   /// Rearranges the rows of the m-by-n matrix X as specified
   /// by the permutation k[0],k[1],...,k[m-1].
   /// If forw is true, a forward permutation is performed:
   ///
   ///        row(k[i]) -> row(i) for each i=0,...,m-1
   /// If forw is false, a backward permutation is performed:
   ///
   ///        row(i) -> row(k[i]) for each i=0,...,m-1
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param forw If true, a forward permutation is performed; if false a backward permutation is performed.
   /// @param m The number of rows of the matrix X.
   /// @param n The number of columns of the matrix X.
   /// @param X Pointer to the complex matrix X. On exit, the rows of X are permuted as specified by k.
   /// @param ldX Column length of the matrix X.  ldX>=m
   /// @param k Pointer to permutation vector of length m.
   /// @ingroup MAT
   
   template<typename real_t>
   int LAPMR(bool forw,int_t m,int_t n,complex<real_t> *X,int_t ldX,int_t *k)
   {
      return LAPMR< complex<real_t> >(forw,m,n,X,ldX,k);
   }
}

#endif
