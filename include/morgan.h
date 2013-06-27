//
//  morgan.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 4/3/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _morgan_h
#define _morgan_h

/// @file morgan.h Generates morgan matrix.

#include <cstdlib>
#include <limits>
#include <cmath>
#include "latl.h"

namespace LATL
{
   /// @brief Generates Morgan matrix of order n and degree k.
   /// @details The Morgan matrix is a lower Hessenberg matrix A where
   ///
   ///           A(i,i)   = i for i=1,...,n
   ///           A(i,i+1) = 1 for i=1,...,n-1
   ///           A(n,1)   = eps^k
   ///           A(i,j)   = 0 otherwise
   /// where eps is the machine epsilon.  This matrix with n=4 and k=2 was reported
   /// by Ron Morgan on the
   /// <A HREF="https://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=13&t=4270">LAPACK forum</A>
   /// to cause the nonsymmetric eigenvalue routine to produce incorrect eigenvectors
   /// unless balancing is disabled.
   /// @tparam T Type of matrix elements.
   /// @return Pointer to matrix of type T containing n-by-n Morgan matrix.
   /// @return nullptr if a parameter is invalid.
   /// @param n Order of matrix.  n>0.
   /// @param k Degree of the epsilon term.
   /// @ingroup MATGEN

   template<typename T> T *Morgan(int_t n,int_t k)
   {
      const T zero(0.0);
      const T one(1.0);
      const T eps=std::numeric_limits<T>::epsilon();
      using std::pow;
      
      if(n<1)
         return nullptr;

      T *M=new T[n*n];
      if(M)
      {
         for(int_t i=0;i<n;i++)
            for(int_t j=0;j<n;j++)
               M[i+j*n]=zero;

         for(int_t i=0;i<n;i++)
            M[i+i*n]=static_cast<T>(i+1);
         
         for(int_t i=0;i<n-1;i++)
            M[i+(i+1)*n]=one;
         
         M[n-1]=pow(eps,k);
      }
      return M;
   }
}
#endif
