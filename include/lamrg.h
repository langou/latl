//
//  lamrg.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/24/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lamrg_h
#define _lamrg_h

/// @file lamrg.h Creates permutation vector for sorting two contiguous lists.

#include "latl.h"

namespace latl
{
   /// @brief Creates permutation vector for sorting two contiguous lists.
   /// 
   /// A permutation list is created which will merge the elements  of A 
   /// (which is composed of two independently sorted sets) into a
   /// single set which is sorted in ascending order.
   /// @tparam real_t Floating point type.
   /// @param n1 Length of the first list stored in a.
   /// @param n2 Length of the second list stored in a.
   /// @param a Pointer to two lists, stored consecutively, stored in either ascending or descending order.
   /// @param s1 True if the first list is stored in ascending order, false if it is stored in descending order.
   /// @param s2 True if the second list is stored in ascending order, false if it is stored in descending order.   
   /// @param index Pointer to permutation list of length n1+n2; on exit this array will contain a permutation
   /// such that if
   ///
   ///        b[i] = a[index[i]]
   /// for i=0,..,n1+n2-1, then b will be sorted in ascending order.
   /// @ingroup VEC
   
   template<typename real_t>
   void lamrg(int_t n1, int_t n2, real_t *a, bool s1, bool s2, int_t *index)
   {
      int_t n1sv=n1;
      int_t n2sv=n2;
      int_t ind1=(s1)?0:n1-1;
      int_t ind2=(s2)?n1:n1+n2-1;
      int_t i=0;
      while((n1sv>0)&&(n2sv>0))
      {
         if(a[ind1]<=a[ind2])
         {
            index[i]=ind1;
            i++;
            if(s1)
               ind1++;
            else
               ind1--;
            n1sv--;
         }
         else
         {
            index[i]=ind2;
            i++;
            if(s2)
               ind2++;
            else
               ind2--;
            n2sv--;
         }
      }
      if(n1sv==0)
      {
         for(int_t j=0;j<n2sv;j++)
         {
            index[i]=ind2;
            i++;
            if(s2)
               ind2++;
            else
               ind2--;
         }
      }
      else
      {
         for(int_t j=0;j<n1sv;j++)
         {
            index[i]=ind1;
            i++;
            if(s1)
               ind1++;
            else
               ind1--;
         }
      }
   }
}

#endif
