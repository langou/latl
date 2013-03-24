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

using namespace LATL;
using namespace std;

#if defined(FLOAT)
typedef float Real;
#elif defined(DOUBLE)
typedef double Real;
#elif defined(LDOUBLE)
typedef long double Real;
#elif defined(REAL)
#include "real.hpp"
typedef mpfr::real<REAL> Real;
#elif defined(MPREAL)
#include "mpreal.h"
typedef mpfr::mpreal Real;
#else
typedef double Real;
#endif
typedef complex<Real> Complex;

// general matrix inverse test

template <typename real_t,typename matrix_t>
real_t general(int nb,bool resi,bool prnt)
{
   char uplo='a';
   int m,n;
   const matrix_t one(1.0);
   const matrix_t zero(0.0);
   real_t residual;
   matrix_t *A=LOAD<matrix_t>(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   matrix_t *B=new matrix_t[n*n];
   matrix_t *C=new matrix_t[n*n];
   matrix_t *D=new matrix_t[n*n];
   int_t *ipiv=new int_t[n];
   LACPY<real_t>(uplo,n,n,A,n,B,n);
   int info=GETRF<real_t>(n,n,A,n,ipiv,nb);
   info=GETRI<real_t>(n,A,n,ipiv,nb);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      PRINT<real_t>(n,n,A,n);
   LASET<real_t>(uplo,n,n,zero,one,C,n);
   LASET<real_t>(uplo,n,n,zero,one,D,n);
   GEMM<real_t>('n','n',n,n,n,one,A,n,B,n,-one,C,n);
   GEMM<real_t>('n','n',n,n,n,one,B,n,A,n,-one,D,n);
   residual=max(LANGE<real_t>('m',n,n,C,n),LANGE<real_t>('m',n,n,D,n));
   real_t condition=LANGE<real_t>('i',n,n,A,n)*LANGE<real_t>('i',n,n,B,n);
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
real_t symmetric(int nb,char uplo,bool resi,bool prnt)
{
   int m,n;
   const matrix_t one(1.0);
   const matrix_t zero(0.0);
   real_t residual;
   matrix_t *A=LOAD<matrix_t>(m,n);
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
   LACPY<real_t>('a',n,n,A,n,B,n);
   int info=SYTRF<real_t>(uplo,n,A,n,ipiv,block,nb);
   info=SYTRI<real_t>(uplo,n,A,n,ipiv,block,nb);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      PRINT<real_t>(n,n,A,n);
   LASET<real_t>('a',n,n,zero,one,C,n);
   LASET<real_t>('a',n,n,zero,one,D,n);
   SYMM<real_t>('l',uplo,n,n,one,A,n,B,n,-one,C,n);
   SYMM<real_t>('r',uplo,n,n,one,A,n,B,n,-one,D,n);
   residual=max(LANGE<real_t>('m',n,n,C,n),LANGE<real_t>('m',n,n,D,n));
   real_t condition=LANGE<real_t>('i',n,n,A,n)*LANGE<real_t>('i',n,n,B,n);
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
real_t hermitian(int nb,char uplo,bool resi,bool prnt)
{
   int m,n;
   const complex<real_t> one(1.0);
   const complex<real_t> zero(0.0);
   real_t residual;
   complex<real_t> *A=LOAD< complex<real_t> >(m,n);
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
   LACPY<real_t>('a',n,n,A,n,B,n);
   int info=HETRF<real_t>(uplo,n,A,n,ipiv,block,nb);
   info=HETRI<real_t>(uplo,n,A,n,ipiv,block);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      PRINT<real_t>(n,n,A,n);
   LASET<real_t>('a',n,n,zero,one,C,n);
   LASET<real_t>('a',n,n,zero,one,D,n);
   HEMM<real_t>('l',uplo,n,n,one,A,n,B,n,-one,C,n);
   HEMM<real_t>('r',uplo,n,n,one,A,n,B,n,-one,D,n);
   residual=max(LANGE<real_t>('m',n,n,C,n),LANGE<real_t>('m',n,n,D,n));
   real_t condition=LANGE<real_t>('i',n,n,A,n)*LANGE<real_t>('i',n,n,B,n);
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
real_t symmetric_positive(int nb,char uplo,bool resi,bool prnt)
{
   int m,n;
   const real_t one(1.0);
   const real_t zero(0.0);
   real_t residual;
   real_t *A=LOAD<real_t>(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   real_t *B=new real_t[n*n];
   real_t *C=new real_t[n*n];
   real_t *D=new real_t[n*n];
   LACPY<real_t>('a',n,n,A,n,B,n);
   int info=POTRF<real_t>(uplo,n,A,n,nb);
   if(info>0)
   {
      cerr << "input matrix is not positive definite" << endl;
      exit(0);
   }
   info=POTRI<real_t>(uplo,n,A,n,nb);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      PRINT<real_t>(n,n,A,n);
   LASET<real_t>('a',n,n,zero,one,C,n);
   LASET<real_t>('a',n,n,zero,one,D,n);
   SYMM<real_t>('l',uplo,n,n,one,A,n,B,n,-one,C,n);
   SYMM<real_t>('r',uplo,n,n,one,A,n,B,n,-one,D,n);
   residual=max(LANGE<real_t>('m',n,n,C,n),LANGE<real_t>('m',n,n,D,n));
   real_t condition=LANGE<real_t>('i',n,n,A,n)*LANGE<real_t>('i',n,n,B,n);
   real_t eps=numeric_limits<real_t>::epsilon();
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return (resi)? residual : residual/(n*condition*eps);
}

// hermitian positive definite matrix inverse test

template <typename real_t>
real_t hermitian_positive(int nb,char uplo,bool resi,bool prnt)
{
   int m,n;
   const complex<real_t> one(1.0);
   const complex<real_t> zero(0.0);
   real_t residual;
   complex<real_t> *A=LOAD< complex<real_t> >(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   complex<real_t> *B=new complex<real_t>[n*n];
   complex<real_t> *C=new complex<real_t>[n*n];
   complex<real_t> *D=new complex<real_t>[n*n];
   LACPY<real_t>('a',n,n,A,n,B,n);
   int info=POTRF<real_t>(uplo,n,A,n,nb);
   if(info>0)
   {
      cerr << "input matrix is not positive definite" << endl;
      exit(0);
   }
   info=POTRI<real_t>(uplo,n,A,n,nb);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      PRINT<real_t>(n,n,A,n);
   LASET<real_t>('a',n,n,zero,one,C,n);
   LASET<real_t>('a',n,n,zero,one,D,n);
   HEMM<real_t>('l',uplo,n,n,one,A,n,B,n,-one,C,n);
   HEMM<real_t>('r',uplo,n,n,one,A,n,B,n,-one,D,n);
   residual=max(LANGE<real_t>('m',n,n,C,n),LANGE<real_t>('m',n,n,D,n));
   real_t condition=LANGE<real_t>('i',n,n,A,n)*LANGE<real_t>('i',n,n,B,n);
   real_t eps=numeric_limits<real_t>::epsilon();
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return (resi)? residual : residual/(n*condition*eps);
}

// triangular matrix inverse test

template <typename real_t,typename matrix_t>
real_t triangular(int nb,char uplo, char diag,bool resi,bool prnt)
{
   int m,n;
   const matrix_t one(1.0);
   const matrix_t zero(0.0);
   real_t residual;
   matrix_t *A=LOAD<matrix_t>(m,n);
   if((A==NULL)||(m!=n))
   {
      cerr << "error reading input matrix" << endl;
      exit(0);
   }
   matrix_t *B=new matrix_t[n*n];
   matrix_t *C=new matrix_t[n*n];
   matrix_t *D=new matrix_t[n*n];
   LACPY<real_t>(uplo,n,n,A,n,B,n);
   int info=TRTRI<real_t>(uplo,diag,n,A,n,nb);
   if(info>0)
   {
      cerr << "input matrix is singular" << endl;
      exit(0);
   }
   if(prnt)
      PRINT<real_t>(uplo,diag,n,A,n);
   LASET<real_t>(uplo,n,n,zero,one,C,n);
   LASET<real_t>(uplo,n,n,zero,one,D,n);
   if(diag=='u')
      for(int i=0;i<n;i++)
         A[i+i*n]=B[i+i*n]=one;
   TTMM<real_t>(uplo,n,one,A,n,B,n,-one,C,n);
   TTMM<real_t>(uplo,n,one,B,n,A,n,-one,D,n);
   residual=max(LANTR<real_t>('m',uplo,'n',n,n,C,n),LANTR<real_t>('m',uplo,'n',n,n,D,n));
   real_t condition=LANTR<real_t>('i',uplo,diag,n,n,A,n)*LANTR<real_t>('i',uplo,diag,n,n,B,n);
   real_t eps=numeric_limits<real_t>::epsilon();
   delete [] D;
   delete [] C;
   delete [] B;
   delete [] A;
   return (resi)? residual : residual/(n*condition*eps);
}

void usage(char *name,int nb)
{
   cerr << "Usage: " << name << " [-g | -t | -s | -h | -p]";
#ifdef MPREAL
   cerr << " [-P <n>]";
#endif
   cerr << " [-c] [-l] [-u [-w] [-r] [-b <nb>]" << endl;
   cerr << "           -g      invert general matrix (default)" << endl;
   cerr << "           -t      invert upper or lower triangular matrix" << endl;
   cerr << "           -s      invert symmetric matrix" << endl;
   cerr << "           -h      invert hermitian matrix" << endl;
   cerr << "           -p      invert positive definite matrix" << endl;
   cerr << "           -c      use complex (default is to use real)" << endl;
   cerr << "           -l      use lower triangular matrix (default is upper triangular)" << endl;
   cerr << "           -u      assume unit triangular matrix" << endl;
   cerr << "           -w      write inverse matrix to standard output" << endl;
   cerr << "           -r      report the residual error instead of the relative error" << endl;
   cerr << "           -b <nb> use block size of nb, otherwise blocksize is set to " << nb << endl;
#ifdef MPREAL
   cerr << "           -P <n>  set mpreal precision to n bits " << endl;
#endif
}

int main(int argc,char **argv)
{

   int arg=1;
   char uplo='u';
   char diag='n';
   char Type='G';
   bool comp=0;
   bool resi=0;
   bool prnt=0;
   int nb=64;
   
   while(arg<argc)
   {
      if(strncmp(argv[arg],"-t",2)==0)
         Type='T';
      else if(strncmp(argv[arg],"-g",2)==0)
         Type='G';
      else if(strncmp(argv[arg],"-s",2)==0)
         Type='S';
      else if(strncmp(argv[arg],"-h",2)==0)
         Type='H';
      else if(strncmp(argv[arg],"-p",2)==0)
         Type='P';
      else if(strncmp(argv[arg],"-l",2)==0)
         uplo='l';
      else if(strncmp(argv[arg],"-c",2)==0)
         comp=1;
      else if(strncmp(argv[arg],"-u",2)==0)
         diag='u';
      else if(strncmp(argv[arg],"-w",2)==0)
         prnt=1;
      else if(strncmp(argv[arg],"-r",2)==0)
         resi=1;
      else if(strncmp(argv[arg],"-b",2)==0)
      {
         arg++;
         if(arg<argc)
            nb=atoi(argv[arg]);
         else
            nb=0;
      }
#ifdef MPREAL
      else if(strncmp(argv[arg],"-P",4)==0)
      {
         int prec=53;
         arg++;
         if(arg<argc)
            prec=atoi(argv[arg]);
         mpfr::mpreal::set_default_prec(prec);
      }
#endif
      else
      {
         usage(argv[0],nb);
         return 1;
      }
      arg++;
   }

   Real error;

   if(Type=='T')
   {
      if(comp)
         error=triangular<Real,Complex>(nb,uplo,diag,resi,prnt);
      else
         error=triangular<Real,Real>(nb,uplo,diag,resi,prnt);
      cerr << error << endl;
      return 0;
   }
   else if(Type=='G')
   {
      if(comp)
         error=general<Real,Complex>(nb,resi,prnt);
      else
         error=general<Real,Real>(nb,resi,prnt);
      cerr << error << endl;
      return 0;
   }
   else if(Type=='S')
   {
      if(comp)
         error=symmetric<Real,Complex>(nb,uplo,resi,prnt);
      else
         error=symmetric<Real,Real>(nb,uplo,resi,prnt);
      cerr << error << endl;
      return 0;
   }
   else if(Type=='H')
   {
      error=hermitian<Real>(nb,uplo,resi,prnt);
      cerr << error << endl;
      return 0;
   }
   else if(Type=='P')
   {
      if(comp)
         error=hermitian_positive<Real>(nb,uplo,resi,prnt);
      else
         error=symmetric_positive<Real>(nb,uplo,resi,prnt);
      cerr << error << endl;
      return 0;
   }

   usage(argv[0],nb);
   return 0;
}
