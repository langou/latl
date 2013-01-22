//
//  load.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 8/30/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _load_h
#define _load_h

/// @file load.h Reads text formatted matrix from standard input.

#include <sstream>
#include <iostream>
#include <queue>
#include <algorithm>

namespace latl
{
   /// @brief Attempts to read a text formatted matrix from standard input.
   /// @return Pointer to matrix read from standard input.
   /// @tparam real_t Floating point type.
   /// @param[out] m Number of rows in matrix.
   /// @param[out] n Number of columns in matrix.
   
   template <typename real_t>
   real_t *load(int &m,int &n)
   {
      using std::string;
      using std::queue;
      using std::getline;
      using std::istringstream;
      using std::max;
      string input;
      queue<string> q;
      n=0;
      while(getline(cin,input))
      {
         q.push(input);
         istringstream s(input);
         real_t z;
         int N=0;
         while(s>>z)
            N++;
         n=max(n,N);
      }
      m=q.size();
      real_t *A=new real_t[m*n];
      int i=0;
      while(!q.empty())
      {
         istringstream ist(q.front());
         real_t x;
         int j=0;
         while(ist>>x)
         {
            A[i+j*m]=x;
            j++;
         }
         i++;
         q.pop();
      }
      return A;
   }
}
#endif
