//
//  matgen.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 4/3/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include <complex>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <cstring>
#include <print.h>
#include <hilbert.h>
#include <morgan.h>
#include <grcar.h>
#include <general.h>
#include <symmetric.h>
#include <hermitian.h>
#include <positive.h>

typedef long double ldouble;

#ifndef REAL
#define REAL double
#endif

typedef REAL Real;
typedef std::complex<Real> Complex;

template<typename T> void output(int m,int n,T *A,char *outfile)
{
   using std::cerr;
   using std::endl;
   using LATL::Print;
   
   if(A)
   {
      if(outfile==NULL)
      {
         Print<T>(m,n,A,m);
      }
      else
      {
         std::ofstream out(outfile);
         if(out)
            Print<T>(m,n,A,m,out);
         else
         {
            cerr << "Cannot open output file " << outfile << "." << endl;
            exit(0);
         }
      }
   }
   else
   {
      cerr << "Cannot create matrix." << endl;
      exit(0);
   }
}

template<typename T> void output(int m,int n,complex<T> *A,char *outfile)
{
   using std::cerr;
   using std::endl;
   using LATL::Print;

   if(A)
   {
      if(outfile==NULL)
      {
         Print<T>(m,n,A,m);
      }
      else
      {
         std::ofstream out(outfile);
         if(out)
            Print<T>(m,n,A,m,out);
         else
         {
            cerr << "Cannot open output file " << outfile << "." << endl;
            exit(0);
         }
      }
   }
   else
   {
      cerr << "Cannot create matrix." << endl;
      exit(0);
   }
}

void usage(char *name)
{
   using std::cerr;
   using std::endl;
   
   cerr << "Usage: " << name << " <options>" << endl << endl;;
   cerr << "        -o <file>  save matrix to <file>, otherwise writes matrix to standard output" << endl;
   cerr << "        -m <m>     sets number of rows (default is m=1)" << endl;
   cerr << "        -n <n>     sets number of columns (default is n=1)" << endl;
   cerr << "        -c         use complex numbers if possible" << endl;
   cerr << "        -hilbert   creates n-by-n Hilbert matrix" << endl;
   cerr << "        -morgan    creates n-by-n Morgan matrix" << endl;
   cerr << "        -grcar     creates n-by-n Grcar matrix with k superdiagonals" << endl;
   cerr << "        -random    creates m-by-n random matrix (default)" << endl;
   cerr << "        -symmetric creates n-by-n random symmetric matrix" << endl;
   cerr << "        -hermitian creates n-by-n random Hermitian matrix" << endl;
   cerr << "        -positive  creates n-by-n random positive definite matrix" << endl;
   cerr << "        -help      display help" << endl;
   exit(0);
}


int main(int argc, char** argv)
{
   using std::cout;
   using std::cerr;
   using std::endl;
   using std::strncmp;
   using std::atoi;
   
   using LATL::Print;
   using LATL::Hilbert;
   using LATL::Morgan;
   using LATL::Grcar;
   using LATL::General;
   using LATL::GeneralComplex;
   using LATL::Symmetric;
   using LATL::SymmetricComplex;
   using LATL::Hermitian;
   using LATL::Positive;
   using LATL::PositiveComplex;
   
   enum MatrixType {hilbert,morgan,grcar,random};
   enum SymmetryType {symmetric,hermitian,positive,none};
   
   int m=1;
   int n=1;
   int k=1;
   int arg=1;
   char *outfile=nullptr;
   bool cmplx=false;
   
   MatrixType matrix=random;
   SymmetryType symmetry=none;
   
   while(arg<argc)
   {
      if(strncmp(argv[arg],"-morgan",7)==0)
      {
         matrix=morgan;
      }
      else if(strncmp(argv[arg],"-hilbert",8)==0)
      {
         matrix=hilbert;
      }
      else if(strncmp(argv[arg],"-grcar",6)==0)
      {
         matrix=grcar;
      }
      else if(strncmp(argv[arg],"-random",7)==0)
      {
         matrix=random;
      }
      else if(strncmp(argv[arg],"-symmetric",4)==0)
      {
         symmetry=symmetric;
      }
      else if(strncmp(argv[arg],"-hermitian",4)==0)
      {
         symmetry=hermitian;
      }
      else if(strncmp(argv[arg],"-positive",4)==0)
      {
         symmetry=positive;
      }
      else if(strncmp(argv[arg],"-help",5)==0)
      {
         usage(argv[0]);
      }
      else if(strncmp(argv[arg],"-c",2)==0)
      {
         cmplx=true;
      }
      else if(strncmp(argv[arg],"-m",2)==0)
      {
         arg++;
         if(arg<argc)
            m=atoi(argv[arg]);
         if(m<1)
            usage(argv[0]);
      }
      else if(strncmp(argv[arg],"-n",2)==0)
      {
         arg++;
         if(arg<argc)
            n=atoi(argv[arg]);
         if(n<1)
            usage(argv[0]);
      }
      else if(strncmp(argv[arg],"-k",2)==0)
      {
         arg++;
         if(arg<argc)
            k=atoi(argv[arg]);
      }
      else if(strncmp(argv[arg],"-o",2)==0)
      {
         arg++;
         if(arg<argc)
            outfile=argv[arg];
      }
      else
      {
         usage(argv[0]);
      }
      arg++;
   }

   if((symmetry==hermitian)&&(!cmplx))
      symmetry=symmetric;

   if(matrix==hilbert)
   {
      Real *A=Hilbert<Real>(n);
      output<Real>(n,n,A,outfile);
   }
   else if(matrix==morgan)
   {
      Real *A=Morgan<Real>(n,k);
      output<Real>(n,n,A,outfile);
   }
   else if(matrix==grcar)
   {
      Real *A=Grcar<Real>(n,k);
      output<Real>(n,n,A,outfile);
   }
   else if(matrix==random)
   {
      if(cmplx)
      {
         if(symmetry==symmetric)
         {
            Complex *A=SymmetricComplex<Real>(n);
            output<Real>(n,n,A,outfile);
         }
         else if(symmetry==hermitian)
         {
            Complex *A=Hermitian<Real>(n);
            output<Real>(n,n,A,outfile);
         }
         else if(symmetry==positive)
         {
            Complex *A=PositiveComplex<Real>(n);
            output<Real>(n,n,A,outfile);
         }
         else
         {
            Complex *A=GeneralComplex<Real>(m,n);
            output<Real>(m,n,A,outfile);
         }
      }
      else
      {
         if(symmetry==symmetric)
         {
            Real *A=Symmetric<Real>(n);
            output<Real>(n,n,A,outfile);
         }
         else if(symmetry==positive)
         {
            Real *A=Positive<Real>(n);
            output<Real>(n,n,A,outfile);
         }
         else
         {
            Real *A=General<Real>(m,n);
            output<Real>(m,n,A,outfile);
         }
      }
   }
   return 0;
}
