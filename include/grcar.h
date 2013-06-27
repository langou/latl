//
//  grcar.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 4/3/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _grcar_h
#define _grcar_h

/// @file grcar.h Generates grcar matrix.

#include <cstdlib>
#include <limits>
#include <cmath>
#include "latl.h"

namespace LATL
{
   /// @brief Generates grcar matrix of order n and k superdiagonals.
   /// @details The grcar matrix is a n-by-n Toeplitz matrix with A with
   ///
   ///           A(i,i)   =  1 for i=1,...,n
   ///           A(i,i-1) = -1 for i=2,...,n
   ///           A(i,i+1) =  1 for i=1,...,n-1
   ///                    .
   ///                    .
   ///                    .
   ///           A(i,i+k) =  1 for i=1,...,n-k
   ///           A(i,j)   =  0 otherwise
   /// Computing eigenvalues of this matrix is very sensitive to the numerical precision used.
   /// Reference: <A HREF="http://math.nist.gov/MatrixMarket/data/NEP/mvmgrc/mvmgrc.html">J. Grcar</A>.
   /// @tparam T Type of matrix elements.
   /// @return Pointer to matrix of type T containing n-by-n grcar matrix.
   /// @return NULL if a parameter is invalid.
   /// @param n Order of matrix.  n>0.
   /// @param k Number of superdiagonals. k>0.
   /// @ingroup MATGEN

   template<typename T> T *Grcar(int_t n,int_t k)
   {
      const T one(1.0);
      const T zero(0.0);
      
      if((n<1)||(k<0))
         return NULL;

      T *M=new T[n*n];
      if(M)
      {
         for(int_t i=0;i<n;i++)
            for(int_t j=0;j<n;j++)
               M[i+j*n]=zero;

         for(int_t i=0;i<n;i++)
            M[i+i*n]=one;
         
         for(int_t i=1;i<n;i++)
            M[i+(i-1)*n]=-one;

         for(int_t j=1;j<=k;j++)
            for(int_t i=0;i<n-j;i++)
               M[i+(i+j)*n]=one;
      }
      return M;
   }
}
#endif
