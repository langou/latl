//
//  matgen.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/24/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//
//  [generates a matrix of specified size and prints to standard output]
//

#include <iostream>
#include <cstdlib>
#include <limits>
#include <cstring>
#include <print.h>
#include <latl.h>
#include <larnv.h>

typedef long double ldouble;

#ifndef REAL
#define REAL double
#endif
typedef REAL real_t;

void usage(char *name)
{
   using std::cerr;
   using std::endl;
   cerr << "Usage: " << name << " [-complex] [-random <dist>] [-m <m>] [-n <n>] [-seed <seed>] [-hilbert] [-symmetric] [-hermitian]" << endl;
   cerr << "        -complex generates complex matrix (default is real)" << endl;
   cerr << "        -m <m> sets number of rows (default is m=1)" << endl;
   cerr << "        -n <n> sets number of columns (default is n=1)" << endl;
   cerr << "        -hilbert creates n-by-n Hilbert matrix" << endl;
   cerr << "        -hermitian creates n-by-n Hermitian matrix" << endl;
   cerr << "        -symmetric creates n-by-n symmetric matrix" << endl;
   cerr << "        -seed <seed> sets random number generator seed (default is 0)" << endl;
   cerr << "        -random <dist> creates M x N random matrix with one of the following distributions: " << endl;
   cerr << "                   1 = uniform on (0,1) (default)" << endl;
   cerr << "                   2 = uniform on (-1,1)" << endl;
   cerr << "                   3 = normal on (0,1)" << endl;
   cerr << "                   4 = uniformly distributed on the disc abs(z) < 1 (complex)" << endl;
   cerr << "                   5 = uniformly distributed on the circle abs(z) = 1 (complex)" << endl;
}

int main(int argc, char** argv)
{
   using std::cout;
   using std::cerr;
   using std::endl;
   using std::numeric_limits;
   using std::strncmp;
   using std::atoi;
   using std::real;
   using std::conj;
   using latl::larnv;
   int m=1;
   int n=1;
   int arg=1;
   uint32_t s=0;
   int dist=1;
   bool use_complex=0;
   bool symmetric=0;
   bool hermitian=0;
   bool hilbert=0;
   
   while(arg<argc)
   {
      if(strncmp(argv[arg],"-m",2)==0)
      {
         arg++;
         m=atoi(argv[arg]);
      }
      else if(strncmp(argv[arg],"-n",2)==0)
      {
         arg++;
         n=atoi(argv[arg]);
      }
      else if(strncmp(argv[arg],"-seed",5)==0)
      {
         arg++;
         s=atoi(argv[arg]);
      }
      else if(strncmp(argv[arg],"-random",7)==0)
      {
         arg++;
         dist=atoi(argv[arg]);
      }
      else if(strncmp(argv[arg],"-complex",8)==0)
      {
         use_complex=1;
      }
      else if(strncmp(argv[arg],"-symmetric",10)==0)
      {
         symmetric=1;
      }
      else if(strncmp(argv[arg],"-hermitian",10)==0)
      {
         hermitian=1;
      }
      else if(strncmp(argv[arg],"-hilbert",8)==0)
      {
         hilbert=1;
      }
      else
      {
         usage(argv[0]);
         return 1;
      }
      arg++;
   }

   if(hilbert)
   {
      hermitian=0;
      symmetric=1;
      use_complex=0;
      dist=0;
   }

   if(hermitian||symmetric)
      m=n;
   
   if(use_complex)
   {
      if(hermitian&&symmetric)
         symmetric=0;
      complex<real_t> *A=new complex<real_t>[m*n];
      if((dist>0)&&(dist<6))
      {
         larnv(dist,m*n,A,s);
         if(symmetric||hermitian)
         {
            for(int i=0;i<n;i++)
            {
               A[i+i*n]=real(A[i+i*n]);
               if(symmetric)
                  for(int j=i+1;j<n;j++)
                     A[i+j*n]=A[j+i*n];
               else
                  for(int j=i+1;j<n;j++)
                     A[i+j*n]=conj(A[j+i*n]);
            }
         }
         latl::print(m,n,A,m);
      }
      else
      {
         usage(argv[0]);
         return 0;
      }
   }
   else
   {
      real_t *A=new real_t[m*n];
      if(hermitian||symmetric)
      {
         hermitian=0;
         symmetric=1;
      }
      if(hilbert)
      {
         for(int j=0;j<n;j++)
            for(int i=0;i<n;i++)
               A[i+j*n]=1.0/(real_t)(i+j+1);
         latl::print(m,n,A,m);
      }
      else if((dist>0)&&(dist<4))
      {
         larnv(dist,m*n,A,s);
         if(symmetric)
         {
            for(int i=0;i<n;i++)
            {
               A[i+i*n]=real(A[i+i*n]);
               for(int j=i+1;j<n;j++)
                  A[i+j*n]=A[j+i*n];
            }
         }
         latl::print(m,n,A,m);
      }
      else
      {
         usage(argv[0]);
         return 0;
      }
   }
   return 0;
}
