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
#include <print.h>
#include <load.h>
#include <ttmm.h>
#include <gemm.h>
#include <symm.h>
#include <hemm.h>
#include <lange.h>
#include <lacpy.h>
#include <laset.h>
#include <lantr.h>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <trtri.h>
#include <getrf.h>
#include <getri.h>
#include <sytrf.h>
#include <sytri.h>
#include <potrf.h>
#include <potri.h>
#include <hetrf.h>
#include <hetri.h>

using namespace latl;
using namespace std;

#ifndef REAL
#define REAL double
#endif

typedef long double ldouble;

// general matrix inverse test

template <typename real_t,typename matrix_t>
real_t general(bool resi,bool prnt)
{
   char uplo='a';
   int m,n;
   const matrix_t one(1.0);
   const matrix_t zero(0.0);
   real_t residual;
   matrix_t *A=load<matrix_t>(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   matrix_t *B=new matrix_t[n*n];
   matrix_t *C=new matrix_t[n*n];
   matrix_t *D=new matrix_t[n*n];
   int_t *ipiv=new int_t[n];
   lacpy<real_t>(uplo,n,n,A,n,B,n);
   int info=getrf<real_t>(n,n,A,n,ipiv);
   info=getri<real_t>(n,A,n,ipiv);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      print<real_t>(n,n,A,n);
   laset<real_t>(uplo,n,n,zero,one,C,n);
   laset<real_t>(uplo,n,n,zero,one,D,n);
   gemm<real_t>('n','n',n,n,n,one,A,n,B,n,-one,C,n);
   gemm<real_t>('n','n',n,n,n,one,B,n,A,n,-one,D,n);
   residual=max(lange<real_t>('m',n,n,C,n),lange<real_t>('m',n,n,D,n));
   real_t condition=lange<real_t>('i',n,n,A,n)*lange<real_t>('i',n,n,B,n);
   real_t eps=numeric_limits<real_t>::epsilon();
   delete [] ipiv;
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return (resi)? residual : residual/(n*condition*eps);
}

// symmetric matrix inverse test

template <typename real_t,typename matrix_t>
real_t symmetric(char uplo,bool resi,bool prnt)
{
   int m,n;
   const matrix_t one(1.0);
   const matrix_t zero(0.0);
   real_t residual;
   matrix_t *A=load<matrix_t>(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   matrix_t *B=new matrix_t[n*n];
   matrix_t *C=new matrix_t[n*n];
   matrix_t *D=new matrix_t[n*n];
   int_t *ipiv=new int_t[n];
   bool *block=new bool[n];
   lacpy<real_t>('a',n,n,A,n,B,n);
   int info=sytf2<real_t>(uplo,n,A,n,ipiv,block);
   info=sytri<real_t>(uplo,n,A,n,ipiv,block);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      print<real_t>(n,n,A,n);
   laset<real_t>('a',n,n,zero,one,C,n);
   laset<real_t>('a',n,n,zero,one,D,n);
   symm<real_t>('l',uplo,n,n,one,A,n,B,n,-one,C,n);
   symm<real_t>('r',uplo,n,n,one,A,n,B,n,-one,D,n);
   residual=max(lange<real_t>('m',n,n,C,n),lange<real_t>('m',n,n,D,n));
   real_t condition=lange<real_t>('i',n,n,A,n)*lange<real_t>('i',n,n,B,n);
   real_t eps=numeric_limits<real_t>::epsilon();
   delete [] block;
   delete [] ipiv;
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return (resi)? residual : residual/(n*condition*eps);
}

// hermitian matrix inverse test

template <typename real_t>
real_t hermitian(char uplo,bool resi,bool prnt)
{
   int m,n;
   const complex<real_t> one(1.0);
   const complex<real_t> zero(0.0);
   real_t residual;
   complex<real_t> *A=load< complex<real_t> >(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   complex<real_t> *B=new complex<real_t>[n*n];
   complex<real_t> *C=new complex<real_t>[n*n];
   complex<real_t> *D=new complex<real_t>[n*n];
   int_t *ipiv=new int_t[n];
   bool *block=new bool[n];
   lacpy<real_t>('a',n,n,A,n,B,n);
   int info=sytf2<real_t>(uplo,n,A,n,ipiv,block);
   info=sytri<real_t>(uplo,n,A,n,ipiv,block);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      print<real_t>(n,n,A,n);
   laset<real_t>('a',n,n,zero,one,C,n);
   laset<real_t>('a',n,n,zero,one,D,n);
   hemm<real_t>('l',uplo,n,n,one,A,n,B,n,-one,C,n);
   hemm<real_t>('r',uplo,n,n,one,A,n,B,n,-one,D,n);
   residual=max(lange<real_t>('m',n,n,C,n),lange<real_t>('m',n,n,D,n));
   real_t condition=lange<real_t>('i',n,n,A,n)*lange<real_t>('i',n,n,B,n);
   real_t eps=numeric_limits<real_t>::epsilon();
   delete [] block;
   delete [] ipiv;
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return (resi)? residual : residual/(n*condition*eps);
}

// symmetric positive definite matrix inverse test

template <typename real_t>
real_t symmetric_positive(char uplo,bool resi,bool prnt)
{
   int m,n;
   const real_t one(1.0);
   const real_t zero(0.0);
   real_t residual;
   real_t *A=load<real_t>(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   real_t *B=new real_t[n*n];
   real_t *C=new real_t[n*n];
   real_t *D=new real_t[n*n];
   lacpy<real_t>('a',n,n,A,n,B,n);
   int info=potrf<real_t>(uplo,n,A,n);
   info=potri<real_t>(uplo,n,A,n);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      print<real_t>(n,n,A,n);
   laset<real_t>('a',n,n,zero,one,C,n);
   laset<real_t>('a',n,n,zero,one,D,n);
   symm<real_t>('l',uplo,n,n,one,A,n,B,n,-one,C,n);
   symm<real_t>('r',uplo,n,n,one,A,n,B,n,-one,D,n);
   residual=max(lange<real_t>('m',n,n,C,n),lange<real_t>('m',n,n,D,n));
   real_t condition=lange<real_t>('i',n,n,A,n)*lange<real_t>('i',n,n,B,n);
   real_t eps=numeric_limits<real_t>::epsilon();
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return (resi)? residual : residual/(n*condition*eps);
}

// hermitian positive definite matrix inverse test

template <typename real_t>
real_t hermitian_positive(char uplo,bool resi,bool prnt)
{
   int m,n;
   const complex<real_t> one(1.0);
   const complex<real_t> zero(0.0);
   real_t residual;
   complex<real_t> *A=load< complex<real_t> >(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   complex<real_t> *B=new complex<real_t>[n*n];
   complex<real_t> *C=new complex<real_t>[n*n];
   complex<real_t> *D=new complex<real_t>[n*n];
   lacpy<real_t>('a',n,n,A,n,B,n);
   int info=potrf<real_t>(uplo,n,A,n);
   info=potri<real_t>(uplo,n,A,n);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      print<real_t>(n,n,A,n);
   laset<real_t>('a',n,n,zero,one,C,n);
   laset<real_t>('a',n,n,zero,one,D,n);
   hemm<real_t>('l',uplo,n,n,one,A,n,B,n,-one,C,n);
   hemm<real_t>('r',uplo,n,n,one,A,n,B,n,-one,D,n);
   residual=max(lange<real_t>('m',n,n,C,n),lange<real_t>('m',n,n,D,n));
   real_t condition=lange<real_t>('i',n,n,A,n)*lange<real_t>('i',n,n,B,n);
   real_t eps=numeric_limits<real_t>::epsilon();
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return (resi)? residual : residual/(n*condition*eps);
}

// triangular matrix inverse test

template <typename real_t,typename matrix_t>
real_t triangular(char uplo, char diag,bool resi,bool prnt)
{
   int m,n;
   const matrix_t one(1.0);
   const matrix_t zero(0.0);
   real_t residual;
   matrix_t *A=load<matrix_t>(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   matrix_t *B=new matrix_t[n*n];
   matrix_t *C=new matrix_t[n*n];
   matrix_t *D=new matrix_t[n*n];
   lacpy<real_t>(uplo,n,n,A,n,B,n);
   int info=trtri<real_t>(uplo,diag,n,A,n);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      print<real_t>(uplo,diag,n,A,n);
   laset<real_t>(uplo,n,n,zero,one,C,n);
   laset<real_t>(uplo,n,n,zero,one,D,n);
   if(diag=='u')
      for(int i=0;i<n;i++)
         A[i+i*n]=B[i+i*n]=one;
   ttmm<real_t>(uplo,n,one,A,n,B,n,-one,C,n);
   ttmm<real_t>(uplo,n,one,B,n,A,n,-one,D,n);
   residual=max(lantr<real_t>('m',uplo,'n',n,n,C,n),lantr<real_t>('m',uplo,'n',n,n,D,n));
   real_t condition=lantr<real_t>('i',uplo,diag,n,n,A,n)*lantr<real_t>('i',uplo,diag,n,n,B,n);
   real_t eps=numeric_limits<real_t>::epsilon();
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return (resi)? residual : residual/(n*condition*eps);
}

void usage(char *name)
{
   cerr << "Usage: " << name << " [-general | -triangular | -symmetric | -hermitian | -positive]";
   cerr << " [-complex] [-lower] [-unit] [-print] [-residual]" << endl;
   cerr << "           -general      invert general matrix (default)" << endl;
   cerr << "           -triangular   invert upper or lower triangular matrix" << endl;
   cerr << "           -symmetric    invert symmetric matrix" << endl;
   cerr << "           -hermitian    invert hermitian matrix" << endl;
   cerr << "           -positive     invert positive definite matrix" << endl;
   cerr << "           -complex      use complex (default is to use real)" << endl;
   cerr << "           -lower        use lower triangular matrix (default is upper triangular)" << endl;
   cerr << "           -unit         assume unit triangular matrix" << endl;
   cerr << "           -print        write inverse matrix to standard output" << endl;
   cerr << "           -residual     report the residual error instead of the relative error" << endl;
}

int main(int argc,char **argv)
{
   typedef REAL Real;
   typedef complex<REAL> Complex;

   int arg=1;
   char uplo='u';
   char diag='n';
   char Type='G';
   bool comp=0;
   bool resi=0;
   bool prnt=0;

   while(arg<argc)
   {
      if(strncmp(argv[arg],"-triangular",2)==0)
         Type='T';
      else if(strncmp(argv[arg],"-general",2)==0)
         Type='G';
      else if(strncmp(argv[arg],"-symmetric",2)==0)
         Type='S';
      else if(strncmp(argv[arg],"-hermitian",2)==0)
         Type='H';
      else if(strncmp(argv[arg],"-positive",3)==0)
         Type='P';
      else if(strncmp(argv[arg],"-lower",2)==0)
         uplo='l';
      else if(strncmp(argv[arg],"-complex",2)==0)
         comp=1;
      else if(strncmp(argv[arg],"-unit",2)==0)
         diag='u';
      else if(strncmp(argv[arg],"-print",3)==0)
         prnt=1;
      else if(strncmp(argv[arg],"-residual",2)==0)
         resi=1;
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
      if(comp)
         error=triangular<Real,Complex>(uplo,diag,resi,prnt);
      else
         error=triangular<Real,Real>(uplo,diag,resi,prnt);
      cerr << error << endl;
      return 0;
   }
   else if(Type=='G')
   {
      if(comp)
         error=general<Real,Complex>(resi,prnt);
      else
         error=general<Real,Real>(resi,prnt);
      cerr << error << endl;
      return 0;
   }
   else if(Type=='S')
   {
      if(comp)
         error=symmetric<Real,Complex>(uplo,resi,prnt);
      else
         error=symmetric<Real,Real>(uplo,resi,prnt);
      cerr << error << endl;
      return 0;
   }
   else if(Type=='H')
   {
      error=hermitian<Real>(uplo,resi,prnt);
      cerr << error << endl;
      return 0;
   }
   else if(Type=='P')
   {
      if(comp)
         error=hermitian_positive<Real>(uplo,resi,prnt);
      else
         error=symmetric_positive<Real>(uplo,resi,prnt);
      cerr << error << endl;
      return 0;
   }

   usage(argv[0]);
   return 0;
}
