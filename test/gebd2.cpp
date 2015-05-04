//
//  gebd2.cpp
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 12/16/14.
//  Copyright (c) 2014 University of Colorado Denver. All rights reserved.
//
//  Bidiagonalizes a matrix and report the relative error
//
// ||A - U * B * V^T||/||A||
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include <print.h>
#include <load.h>
#include "timer.h"

#include <latl.h>
#include <gebd2.h>
#include <gemm.h>
#include <lacpy.h>
#include <lange.h>
#include <laset.h>
#include <org2r.h>
#include <orgl2.h>
#include <ormr2.h>

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

template <typename real_t, typename matrix_t>
void get_left_orthogonal (int m, int n, matrix_t *A, real_t *tauq, matrix_t *Q)
{
  int_t info;
  LASET<real_t>('a', m, m, 0, 1, Q, m);
  if(n<=m)
  {
    // A is upper bidiagonal
    LACPY<real_t>('a', m, n, A, m, Q, m);
    info=ORG2R<real_t>(m, m, n, Q, m, tauq);
  }
  else
  {
    // A is lower bidiagonal
    LACPY<real_t>('a', m-1, m-1, &A[1], m, &Q[1+m], m);
    info=ORG2R<real_t>(m-1, m-1, m-1, &Q[1+m], m, tauq);
  }
  if (info!=0)
  {
    cerr <<"Construction of left orthogonal matrix failed. INFO=" <<info
	 <<endl;
    exit(0);
  }
}

template <typename real_t, typename matrix_t>
void get_right_orthogonal (int m, int n, matrix_t *A, real_t *taup, matrix_t *P)
{
  int_t info;

  LASET<real_t>('a', n, n, 0, 1, P, n);
  if (n<=m)
  {
    // A is upper bidiagonal
    LACPY<real_t>('a', n-1, n-1, &A[m], m, &P[1+n], n);
    info=ORGL2<real_t>(n-1, n-1, n-1, &P[1+n], n, taup);
  }
  else
  {
    // A is lower bidiagonal
    LACPY<real_t>('a', m, n, A, m, P, n);
    info=ORGL2<real_t>(n, n, m, P, n, taup);
  }
  if (info!=0)
  {
    cerr <<"Construction of right orthogonal matrix failed. INFO=" <<info
  	 <<endl;
    exit(0);
  }
}

template <typename real_t, typename matrix_t>
void compute_orthogonal_residues (int m, int n, matrix_t *A,
				  real_t *tauq, real_t *taup, bool prnt)
{
  const int d=max(m,n);
  const real_t one(1.0);
  real_t res1, res2;
  matrix_t *R=new matrix_t[d*d];
  
  matrix_t *Q=new matrix_t[m*m];
  get_left_orthogonal<real_t, matrix_t> (m,n,A,tauq,Q);
  if (prnt)
  {
    cout <<"Left orthogonal factor Q:" <<endl;
    Print<real_t>(m,m,Q,m);
  }
  LASET<real_t>('a', m, m, 0, 1, R, m);
  GEMM<real_t>('n', 't', m, m, m, one, Q, m, Q, m, -one, R, m);
  res1 = LANGE<real_t>('m', m, m, R, m);
  LASET<real_t>('a', m, m, 0, 1, R, m);
  GEMM<real_t>('t', 'n', m, m, m, one, Q, m, Q, m, -one, R, m);
  res2=LANGE<real_t>('m', m, m, R, m);
  cout <<"Residue for left orthogonal factor: " <<max(res1, res2) <<endl;

  matrix_t *P=new matrix_t[n*n];
  get_right_orthogonal<real_t, matrix_t> (m,n,A,taup,P);
  if (prnt)
  {
    cout <<"Right orthogonal factor P:" <<endl;
    Print<real_t>(n,n,P,n);
  }
  LASET<real_t>('a', n, n, 0, 1, R, n);
  GEMM<real_t>('n', 't', n, n, n, one, P, n, P, n, -one, R, n);
  res1=LANGE<real_t>('m', n, n, R, n);
  LASET<real_t>('a', n, n, 0, 1, R, n);
  GEMM<real_t>('t', 'n', n, n, n, one, P, n, P, n, -one, R, n);
  res2=LANGE<real_t>('m', n, n, R, n);
  cout <<"Residue for right orthogonal factor: " <<max(res1,res2) <<endl;

  delete [] R;
  delete [] P;
  delete [] Q;
}

template <typename matrix_t>
void check_bidiagonal (int_t m, int_t n, matrix_t *A,
		       matrix_t *D, matrix_t *E)
{
  const int k=min(m,n);
  int shift;
  if (m>=n)
    shift=m; // vector E is super diagonal
  else
    shift=1; // vector E is subdiagonal

  for(int i=0; i<k; i++)
  {
    if (D[i] != A[i+i*m])
    {
      cerr <<"Vector D does not match diagonal elements.";
      cerr <<"A(" <<i <<"," <<i <<")=" <<A[i+i*m] <<endl;
      cerr <<"D(" <<i <<")=" <<D[i] <<endl;
      exit(0);
    }
  }
  for(int i=0; i<k-1; i++)
  {
    if (E[i] != A[i+i*m+shift])
    {
      cerr <<"Vector E does not match subdiagonal elements.";
      cerr <<"A(" <<i+1 <<"," <<i <<")=" <<A[i+i*m+shift] <<endl;
      cerr <<"E(" <<i <<")=" <<E[i] <<endl;
      exit(0);
    }
  }
}

template <typename real_t, typename matrix_t>
void copy_bidiagonal(int m, int n, matrix_t *A, matrix_t *B)
{
  for (int i=0; i<m; i++)
    for (int j=0; j<n; j++)
      if (i==j
	  || (m<n && i==j+1)
	  || (m>=n && j==i+1))
	B[i+j*m]=A[i+j*m];
      else
	B[i+j*m]=0;
}

template <typename real_t, typename matrix_t>
int test_gebd2(bool prnt, bool res, bool time)
{
  int m, n;
  const real_t zero(0.0);
  const real_t one(1.0);

  matrix_t *A=Load<matrix_t>(m,n);
  if (A==NULL)
  {
    cerr <<"error reading input matrix" << endl;
    exit(0);
  }
  if (m==0 || n==0)
    // empty matrix
    return -1;

  const int k=min(m,n);
  matrix_t *D=new matrix_t[k];
  matrix_t *E=new matrix_t[k-1];
  real_t *tauq=new real_t[k];
  real_t *taup=new real_t[k];

  real_t norm=LANGE<real_t>('m', m, n, A, m);
  matrix_t *AA=new matrix_t[m*n];
  LACPY<real_t> ('a',m,n,A,m,AA,m);

  Timer timeit;
  timeit.Start();
  int info=GEBD2<real_t>(m,n,A,m,D,E,tauq,taup);
  timeit.Stop();
  if (info!=0)
  {
    cout <<"GEBD2 failed. INFO=" <<info <<endl;
    exit(0);
  }

  if (prnt)
  {
    cout <<"output A:" <<endl;
    Print<real_t>(m,n,A,m);
    cout <<"tauq:" <<endl;
    Print<real_t>(1,k,tauq,1);
    cout <<"taup:" <<endl;
    Print<real_t>(1,k,taup,1);
  }

  check_bidiagonal<matrix_t> (m, n, A, D, E);
  if (res)
    compute_orthogonal_residues<real_t,matrix_t> (m,n,A,tauq,taup,prnt);

  matrix_t *B=new matrix_t[m*n];
  copy_bidiagonal<real_t,matrix_t>(m,n,A,B);
  matrix_t *Q=new matrix_t[m*m];
  get_left_orthogonal<real_t, matrix_t> (m,n,A,tauq,Q);
  matrix_t *P=new matrix_t[n*n];
  get_right_orthogonal<real_t, matrix_t> (m,n,A,taup,P);

  GEMM<real_t>('n','n',m,n,m,one,Q,m,B,m,zero,A,m);
  GEMM<real_t>('n','n',m,n,n,one,A,m,P,n,-one,AA,m);

  real_t error=LANGE<real_t>('m',m,n,AA,m);
  cout <<"Error: " <<error <<endl;
  cout <<"Relative error: " <<error/norm <<endl;

  if (time)
    cout <<"Elapsed time: " <<timeit.Time() <<" us" <<endl;
  
  delete [] P;
  delete [] Q;
  delete [] B;
  delete [] AA;
  delete [] taup;
  delete [] tauq;
  delete [] E;
  delete [] D;
  delete [] A;
  return 0;
}

void usage(char *name)
{
  cerr << "Usage: " <<name << " [-w] [-r]";
#ifdef MPREAL
  cerr << " [-P <n>]";
#endif
  cerr << endl;
  cerr << "           -w      write matrices to standard output" << endl;
  cerr << "           -r      write residues to standard output" << endl;
  cerr << "           -t      write elapsed time to standard output" << endl;
#ifdef MPREAL
  cerr << "           -P <n>  set mpreal precision to n bits " << endl;
#endif
}

int main(int argc,char **argv)
{
	using std::strncmp;
	
  int arg=1;
  bool prnt=0;
  bool res=0;
  bool time=0;

  while(arg<argc)
  {
    if(strncmp(argv[arg],"-w",2)==0)
      prnt=1;
    else if(strncmp(argv[arg],"-r",2)==0)
      res=1;
    else if(strncmp(argv[arg],"-t",2)==0)
      time=1;
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
      usage(argv[0]);
      return 1;
    }
    arg++;
  }

  if (test_gebd2<Real,Real>(prnt, res, time) !=0)
    usage(argv[0]);

  return 0;
}
