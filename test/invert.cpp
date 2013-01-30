//
//  invert.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/25/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//
//  Computes inverse of a matrix and reports its residual error:
//
//  max{ || A * A^(-1) - I ||,  || A^(-1) * A - I || } / (cond(A) * epsilon * n)
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
#include <cstring>

using namespace latl;
using namespace std;

#ifndef REAL
#define REAL double
#endif

typedef long double ldouble;

// triangular matrix inverse test

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
   lacpy<real_t>(uplo,n,n,A,n,B,n);
   int info=trti2<real_t>(uplo,diag,n,A,n);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   real_t cond=lantr<real_t>('i',uplo,diag,n,n,A,n)*lantr<real_t>('i',uplo,diag,n,n,B,n);
   laset<real_t>(uplo,n,n,zero,one,C,n);
   laset<real_t>(uplo,n,n,zero,one,D,n);
   if(diag=='u')
      for(int i=0;i<n;i++)
         A[i+i*n]=B[i+i*n]=one;
   ttmm<real_t>(uplo,n,one,A,n,B,n,-one,C,n);
   ttmm<real_t>(uplo,n,one,B,n,A,n,-one,D,n);
   residual=max(lantr<real_t>('m',uplo,'n',n,n,C,n),lantr<real_t>('m',uplo,'n',n,n,D,n));
   real_t eps=numeric_limits<real_t>::epsilon();
//   cout << "condition = " << cond << endl;
//   cout << "epsilon   = " << eps << endl;
//   cout << "residual  = " << residual << endl;
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return residual/(n*cond*eps);
}

void usage(char *name)
{
   cerr << "Usage: " << name << " -t <uplo> [-c] [-u]" << endl;
   cerr << endl;
   cerr << "           -c : use complex (default is to use real)" << endl;
   cerr << "           -t <uplo> : invert triangular matrix where uplo is either U or L" << endl;
   cerr << "           -u : assume unit triangular matrix" << endl;
}

int main(int argc,char **argv)
{
   typedef REAL Real;
   typedef complex<REAL> Complex;

   int arg=1;
   char uplo;
   char diag='n';
   char Type='_';
   bool UseComplex=0;

   while(arg<argc)
   {
      if(strncmp(argv[arg],"-t",2)==0)
      {
         Type='T';
         arg++;
         uplo=tolower(*argv[arg]);
         if((uplo!='l')&&(uplo!='u'))
         {
            cerr << "uplo must be either U for upper or L for lower" << endl;
            return 1;
         }
      }
      else if(strncmp(argv[arg],"-c",2)==0)
      {
         UseComplex=1;
      }
      else if(strncmp(argv[arg],"-u",2)==0)
      {
         diag='u';
      }
      else
      {
         usage(argv[0]);
         return 1;
      }
      arg++;
   }

   Real error;

   if(Type=='T')
   {
      if(UseComplex)
         error=triangular<Real,Complex>(uplo,diag);
      else
         error=triangular<Real,Real>(uplo,diag);
      cout << error << endl;
      return 0;
   }

   usage(argv[0]);
   return 0;
}
