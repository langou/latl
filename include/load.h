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
#include <cstddef>

namespace LATL
{
   /// @brief Attempts to read a rectangular matrix from standard input.
   /// @return Pointer to rectangular matrix read from standard input.
   /// @return NULL if unable to read matrix.
   /// @tparam real_t Floating point type.
   /// @param[out] m Number of rows in matrix.
   /// @param[out] n Number of columns in matrix.
   /// @param[in] infile Input stream from which to read matrix.  Defaults to standard input.
   
   template <typename real_t>
   real_t *load(int &m,int &n,std::istream &infile=std::cin)
   {
      using std::string;
      using std::queue;
      using std::getline;
      using std::istringstream;
      using std::max;

      string input;
      queue<string> q;
      n=0;
      int count=0;
      while(getline(infile,input))
      {
         q.push(input);
         istringstream s(input);
         real_t z;
         int N=0;
         while(s>>z)
            N++;
         count+=N;
         n=max(n,N);
      }
      m=q.size();
      if(m*n!=count)
      {
         m=n=0;
         return NULL;
      }
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

   /// @brief Attempts to read a packed triangular matrix from standard input.
   /// @return Pointer to packed triangular matrix read from standard input.
   /// @return NULL if unable to read matrix.
   /// @tparam real_t Floating point type.
   /// @param[out] uplo Determines whether the matrix is upper or lower triangular.
   ///
   ///             uplo == 'U' : matrix is upper triangular
   ///             uplo == 'L' : matrix is lower triangular
   /// @param[out] n Order of packed triangular matrix.
   /// @param[in] infile Input stream from which to read matrix.  Defaults to standard input.
   
   template <typename real_t>
   real_t *load(char &uplo,int &n,std::istream &infile=std::cin)
   {
      using std::string;
      using std::queue;
      using std::getline;
      using std::istringstream;
      using std::max;
      using std::cin;
      string input;
      queue<string> q;
      n=0;
      int m=0;
      int count=0;
      bool first=1;
      while(getline(infile,input))
      {
         q.push(input);
         istringstream s(input);
         real_t z;
         int N=0;
         while(s>>z)
            N++;
         count+=N;
         n=max(n,N);
         if(first)
         {
            m=n;
            first=0;
         }
      }
      if(n*(n+1)/2!=count)
      {
         n=0;
         return NULL;
      }
      bool upper=(m==n)?1:0;
      uplo=(upper)?'U':'L';
      real_t *A=new real_t[(n*(n+1))/2];
      for(int i=0;i<n;i++)
      {
         istringstream ist(q.front());
         if(upper)
         {
            for(int j=i;j<n;j++)
               ist >> A[i+j*(j+1)/2];
         }
         else
         {
            for(int j=0;j<=i;j++)
            {
               ist >> A[i+j*(2*n-1-j)/2];
            }
         }
         q.pop();
      }
      
      return A;
   }
}
#endif
