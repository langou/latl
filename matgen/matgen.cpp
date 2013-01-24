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
   using latl::larnv;
   int m=1;
   int n=1;
   int arg=1;
   int prec=6;
   uint32_t s=0;
   int dist=1;
   bool comp=0;
   
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
      else
      {
         cerr << "Usage: " << argv[0] << " [-c] [-r <dist>] [-m <M>] [-n <N>] [-p <prec>] [-s <seed>]" << endl;
         cerr << "        -r <dist> random number distribution (default = 1)" << endl;
         cerr << "                    real: " << endl;
         cerr << "                          1 = uniform on (0,1)" << endl;
         cerr << "                          2 = uniform on (-1,1)" << endl;
         cerr << "                          3 = normal on (0,1)" << endl;
         cerr << "                    complex:" << endl;
         cerr << "                          1 = real and imaginary uniform on (0,1)" << endl;
         cerr << "                          2 = real and imaginary uniform on (-1,1)" << endl;
         cerr << "                          3 = real and imaginary normal on (0,1)" << endl;
         cerr << "                          4 = uniformly distributed on the disc abs(z) < 1" << endl;
         cerr << "                          5 = uniformly distributed on the circle abs(z) = 1" << endl;
         cerr << "        -c generate complex matrix (default is real)" << endl;
         cerr << "        -p <prec> set precision to <prec> (default is 6)" << endl;
         cerr << "        -m <M> set number of rows to <M> (default is 1)" << endl;
         cerr << "        -n <N> set number of columns to <N> (default is 1)" << endl;
         cerr << "        -s <seed> set random number generator seed to <seed> (default is 0)" << endl;
         return 1;
      }
      arg++;
   }
   
   if(comp)
   {
      if((dist<1)||(dist>5))
      {
         cerr << "invalid distribution" << endl;
         return 1;
      }
      complex<real_t> *A=new complex<real_t>[m*n];
      larnv(dist,m*n,A,s);
      latl::print(m,n,A,m,prec);
   }
   else
   {
      if((dist<1)||(dist>3))
      {
         cerr << "invalid distribution" << endl;
         return 1;         
      }
      real_t *A=new real_t[m*n];
      larnv(dist,m*n,A,s);
      latl::print(m,n,A,m,prec);
   }
   return 0;
}
