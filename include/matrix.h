//
//  matrix.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 4/30/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _matrix_h
#define _matrix_h

/// @file matrix.h Matrix class.

#include <sstream>
#include <iostream>
#include <queue>
#include <algorithm>
#include <cstddef>

namespace LATL
{
   template <typename T> class Matrix
   {
   public:
      const std::size_t row;
      const std::size_t col;

      Matrix(size_t m,size_t n) : row(m),col(n),len(m),allocated(1) { ptr=new T[m*n]; }

      Matrix(std::size_t m,std::size_t n, T* p) : row(m),col(n),ptr(p),len(m),allocated(1) { }

      Matrix(std::size_t m,std::size_t n,std::size_t l,T* p) : row(m),col(n),ptr(p),len(l),allocated(0) { }

      Matrix(Matrix &A,std::ptrdiff_t i1,std::ptrdiff_t i2,std::ptrdiff_t j1,std::ptrdiff_t j2) :
      row(i2-i1+1),col(j2-j1+1),ptr(A.ptr),len(A.len),allocated(0) { }

      ~Matrix() { if(allocated) delete [] ptr; }

      inline T& operator()(std::ptrdiff_t i,std::ptrdiff_t j) { return(ptr[--i+(--j)*len]); }

   private:
      T *ptr;
      const std::size_t len;
      const bool allocated;
   };
}
#endif
