//
//  LATL.cpp
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

typedef long double ldouble;

#ifndef REAL
#define REAL double
#endif

typedef REAL Real;


template<typename T> void output(int m,int n,T *A,char *outfile)
{
   using std::cerr;
   using std::endl;
   using LATL::Print;
   
   if(A)
   {
      if(outfile==NULL)
      {
         Print<T>(n,n,A,n);
      }
      else
      {
         std::ofstream out(outfile);
         if(out)
            Print<T>(m,n,A,n,out);
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
   
   cerr << "Usage: " << name << " [-o <file>] [-m <m>] [-n <n>] [-k <k>] [-hilbert] [-morgan] [-grcar]";
   cerr << endl;
   cerr << "        -o <file> save matrix to <file>, otherwise writes matrix to standard output" << endl;
   cerr << "        -m <m>    sets number of rows (default is m=1)" << endl;
   cerr << "        -n <n>    sets number of columns (default is n=1)" << endl;
   cerr << "        -hilbert  creates n-by-n Hilbert matrix" << endl;
   cerr << "        -morgan   creates n-by-n Morgan matrix" << endl;
   cerr << "        -grcar    creates n-by-n Grcar matrix with k superdiagonals" << endl;
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
   
   enum MatrixType {hilbert,morgan,grcar,none};

   int m=1;
   int n=1;
   int k=1;
   int arg=1;
   char *outfile=NULL;

   MatrixType matrix=none;
   
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
      else if(strncmp(argv[arg],"-grcar",8)==0)
      {
         matrix=grcar;
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
   else
   {
      usage(argv[0]);
   }
   return 0;
}
