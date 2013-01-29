//
//  invert.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/25/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//
//  Computes inverse of a matrix and its residual
//

#include <iostream>
#include <trtri.h>
#include <print.h>
#include <load.h>
#include <ttmm.h>
#include <lacpy.h>
#include <laset.h>
#include <lantr.h>
#include <limits>
#include <cstdlib>

using namespace latl;
using namespace std;

#ifndef REAL
#define REAL double
#endif

typedef long double ldouble;

// triangular matrix

template <typename real_t,typename matrix_t> real_t triangular(char uplo, char diag)
{
   int m,n;
   const matrix_t one(1.0);
   const matrix_t zero(0.0);
   real_t residual;
   matrix_t *A=load<matrix_t>(m,n);
   if(A==NULL)
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   if(m!=n)
   {
      cerr << "input matrix is not square" << endl;
      exit(0);
   }
   matrix_t *B=new matrix_t[n*n];
   matrix_t *C=new matrix_t[n*n];
   matrix_t *D=new matrix_t[n*n];
   lacpy<real_t>('A',n,n,A,n,B,n);
   if(diag=='u')
      for(int i=0;i<n;i++)
         B[i+i*n]=one;
   int info=trti2<real_t>(uplo,diag,n,A,n);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   else
   {
      real_t cond=lantr<real_t>('i',uplo,diag,n,n,A,n)*lantr<real_t>('i',uplo,diag,n,n,B,n);
      laset<real_t>(uplo,n,n,zero,one,C,n);
      laset<real_t>(uplo,n,n,zero,one,D,n);
      ttmm<real_t>(uplo,n,one,A,n,B,n,-one,C,n);
      ttmm<real_t>(uplo,n,one,B,n,A,n,-one,D,n);
      residual=max(lantr<real_t>('m',uplo,diag,n,n,C,n),lantr<real_t>('m',uplo,diag,n,n,D,n));
      real_t eps=numeric_limits<real_t>::epsilon();
      cout << "condition = " << cond << endl;
      cout << "epsilon   = " << eps << endl;
      cout << "residual  = " << residual << endl;
   }
   delete [] D;
   delete [] C;
   delete [] B;
   return residual;
}

int main(int argc,char **argv)
{
   typedef REAL real_t;
   char uplo='L';
   char diag='N';
//   triangular<REAL,REAL>(uplo,diag);
   triangular<REAL,complex<REAL> >(uplo,diag);
   return 0;
}
