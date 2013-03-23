//
//  lasrt.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/21/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lasrt_h
#define _lasrt_h

/// @file lasrt.h Sorts a vector in increasing or decreasing order.

#include <algorithm>
#include <cctype>
#include "latl.h"

template<typename real_t> 
static bool greater_than(real_t a,real_t b) { return(a>b); }

namespace LATL
{
   /// @brief Sorts a vector in increasing or descreasing order.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param direction Specifies ordering of sort.
   ///
   ///        'I' or 'i' for increasing order
   ///        'D' or 'd' for decreasing order
   /// @param n Length of vector d.  n>=0
   /// @param d Pointer to real vector of length n.
   
   template<typename real_t>
   int lasrt(char direction,int_t n,real_t *d)
   {
      direction=std::toupper(direction);
      if((direction!='I')&&(direction!='D'))
         return -1;
      else if(n<0)
         return -2;
      if(n>1)
      {
         if(direction=='I')
            std::sort<real_t*>(d,d+n);
         else
            std::sort<real_t*>(d,d+n,greater_than<real_t>);
      }
      return 0;
   }
}

#endif
