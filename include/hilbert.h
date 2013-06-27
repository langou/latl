//
//  hilbert.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 4/3/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _hilbert_h
#define _hilbert_h

/// @file hilbert.h Generates Hilbert matrix.

#include <cstdlib>
#include "latl.h"

namespace LATL
{
   /// @brief Generates Hilbert matrix of order n.
   /// @details The Hilbert matrix H is defined as
   ///
   ///           H(i,j) = 1/(i+j-1); i,j=1,...,n.
   /// @tparam T Type of matrix elements.
   /// @return Pointer to matrix of type T containing n-by-n Hilbert matrix.
   /// @return nullptr if parameter invalid.
   /// @param n Order of matrix.  n>0.
   /// @ingroup MATGEN

   template<typename T> T *Hilbert(int_t n)
   {
      if(n<1)
         return nullptr;

      T *M=new T[n*n];
      if(M)
      {
         for(int_t i=0;i<n;i++)
            for(int_t j=0;j<=i;j++)
               M[i+j*n]=M[j+i*n]=static_cast<T>(1)/static_cast<T>(i+j+1);
      }
      return M;
   }
}
#endif
