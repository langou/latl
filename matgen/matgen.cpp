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

#ifndef REAL
#define REAL double
#endif
typedef REAL real_t;

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
   int prec=6;
   uint32_t s=0;
   int dist=1;
   bool comp=0;
   bool sym=0;
   bool her=0;
   
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
      else if(strncmp(argv[arg],"-p",2)==0)
      {
         arg++;
         prec=atoi(argv[arg]);
      }
      else if(strncmp(argv[arg],"-s",2)==0)
      {
         arg++;
         s=atoi(argv[arg]);
      }
      else if(strncmp(argv[arg],"-r",2)==0)
      {
         arg++;
         dist=atoi(argv[arg]);
      }
      else if(strncmp(argv[arg],"-c",2)==0)
      {
         comp=1;
      }
      else if(strncmp(argv[arg],"-S",2)==0)
      {
         sym=1;
      }
      else if(strncmp(argv[arg],"-H",2)==0)
      {
         her=1;
      }
      else
      {
         cerr << "Usage: " << argv[0] << " [-c] [-r <dist>] [-m <M>] [-n <N>] [-p <prec>] [-s <seed>]" << endl;
         cerr << "        -c generate complex matrix (default is real)" << endl;
         cerr << "        -p <prec> set precision to <prec> (default is 6)" << endl;
         cerr << "        -m <M> set number of rows to <M> (default is 1)" << endl;
         cerr << "        -n <N> set number of columns to <N> (default is 1)" << endl;
         cerr << "        -H create Hermitian matrix" << endl;
         cerr << "        -S create symmetric matrix" << endl;
         cerr << "        -s <seed> set random number generator seed to <seed> (default is 0)" << endl;
         cerr << "        -r <dist> random number distribution (default = 1)" << endl;
         cerr << "                  real:    1 = uniform on (0,1)" << endl;
         cerr << "                           2 = uniform on (-1,1)" << endl;
         cerr << "                           3 = normal on (0,1)" << endl;
         cerr << "                  complex: 1 = real and imaginary uniform on (0,1)" << endl;
         cerr << "                           2 = real and imaginary uniform on (-1,1)" << endl;
         cerr << "                           3 = real and imaginary normal on (0,1)" << endl;
         cerr << "                           4 = uniformly distributed on the disc abs(z) < 1" << endl;
         cerr << "                           5 = uniformly distributed on the circle abs(z) = 1" << endl;
         return 1;
      }
      arg++;
   }

   if(her||sym)
      m=n;
   if(comp)
   {
      if((dist<1)||(dist>5))
      {
         cerr << "invalid distribution" << endl;
         return 1;
      }
      if(her&&sym)
         sym=0;
      complex<real_t> *A=new complex<real_t>[m*n];
      larnv(dist,m*n,A,s);
      if(sym||her)
      {
         for(int i=0;i<n;i++)
         {
            A[i+i*n]=real(A[i+i*n]);
            if(sym)
               for(int j=i+1;j<n;j++)
                  A[i+j*n]=A[j+i*n];
            else
               for(int j=i+1;j<n;j++)
                  A[i+j*n]=conj(A[j+i*n]);
         }
      }
      latl::print(m,n,A,m,prec);
   }
   else
   {
      if(her||sym)
      {
         her=0;
         sym=1;
      }
      if((dist<1)||(dist>3))
      {
         cerr << "invalid distribution" << endl;
         return 1;         
      }
      real_t *A=new real_t[m*n];
      larnv(dist,m*n,A,s);
      if(sym)
      {
         for(int i=0;i<n;i++)
         {
            A[i+i*n]=real(A[i+i*n]);
            for(int j=i+1;j<n;j++)
               A[i+j*n]=A[j+i*n];
         }
      }
      latl::print(m,n,A,m,prec);
   }
   return 0;
}
