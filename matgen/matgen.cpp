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
#include <fstream>
#include <cstdlib>
#include <cmath>
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
   
   cerr << "Usage: " << name << " [-o <file>] [-c] [-r <dist>] [-m <m>] [-n <n>] [-S <seed>] [-H] [-s] [-h] [-p]";
   cerr << endl;
   cerr << "        -o <file> save matrix to <file>, otherwise writes matrix to standard output" << endl;
   cerr << "        -c        generates complex matrix (default is real)" << endl;
   cerr << "        -m <m>    sets number of rows (default is m=1)" << endl;
   cerr << "        -n <n>    sets number of columns (default is n=1)" << endl;
   cerr << "        -H        creates n-by-n Hilbert matrix" << endl;
   cerr << "        -h        creates n-by-n hermitian matrix" << endl;
   cerr << "        -s        creates n-by-n symmetric matrix" << endl;
   cerr << "        -p        creates n-by-n positive definite matrix by adding n to the diagonal" << endl;
   cerr << "        -S <seed> sets random number generator seed (default is 0)" << endl;
   cerr << "        -e        creates M x N random matrix with entries exp(x), where x is determined below" << endl;
   cerr << "        -r <dist> creates M x N random matrix with one of the following distributions: " << endl;
   cerr << "                    1 = uniform on (0,1) (default)" << endl;
   cerr << "                    2 = uniform on (-1,1)" << endl;
   cerr << "                    3 = normal on (0,1)" << endl;
   cerr << "                    4 = uniformly distributed on the disc abs(z) < 1 (complex)" << endl;
   cerr << "                    5 = uniformly distributed on the circle abs(z) = 1 (complex)" << endl;
   exit(0);
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
   using std::exp;
   using std::log;
   using std::log2;
   using LATL::LARNV;
   int m=1;
   int n=1;
   int arg=1;
   uint32_t s=0;
   int dist=1;
   bool use_complex=0;
   bool symmetric=0;
   bool hermitian=0;
   bool hilbert=0;
   bool positive=0;
   char *outfile=NULL;
   
   while(arg<argc)
   {
      if(strncmp(argv[arg],"-m",2)==0)
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
      else if(strncmp(argv[arg],"-S",2)==0)
      {
         arg++;
         if(arg<argc)
            s=atoi(argv[arg]);
         if(s<0)
            usage(argv[0]);
      }
      else if(strncmp(argv[arg],"-r",2)==0)
      {
         arg++;
         if(arg<argc)
            dist=atoi(argv[arg]);
         if(dist<0)
            usage(argv[0]);
      }
      else if(strncmp(argv[arg],"-c",2)==0)
      {
         use_complex=1;
      }
      else if(strncmp(argv[arg],"-s",2)==0)
      {
         symmetric=1;
      }
      else if(strncmp(argv[arg],"-h",2)==0)
      {
         hermitian=1;
      }
      else if(strncmp(argv[arg],"-H",2)==0)
      {
         hilbert=1;
      }
      else if(strncmp(argv[arg],"-p",2)==0)
      {
         positive=1;
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
   
   if(hilbert)
   {
      hermitian=0;
      symmetric=1;
      use_complex=0;
      dist=0;
   }

   if(positive)
   {
      if(use_complex)
         hermitian=1;
      else
         symmetric=1;
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
         LARNV(dist,m*n,A,s);
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
         if(positive)
         {
            for(int i=0;i<n;i++)
               A[i+n*i]+=(real_t)n;
         }
      }
      if(outfile==NULL)
      {
         LATL::Print<real_t>(m,n,A,m);
      }
      else
      {
         std::ofstream out(outfile);
         if(out)
            LATL::Print<real_t>(m,n,A,m,out);
         else
         {
            cerr << "Cannot open output file." << endl;
            return 1;
         }
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
      }
      else if((dist>0)&&(dist<4))
      {
         LARNV(dist,m*n,A,s);
         if(symmetric)
         {
            for(int i=0;i<n;i++)
            {
               A[i+i*n]=real(A[i+i*n]);
               for(int j=i+1;j<n;j++)
                  A[i+j*n]=A[j+i*n];
            }
         }
         if(positive)
         {
            for(int i=0;i<n;i++)
               A[i+n*i]+=(real_t)n;
         }
      }
      if(outfile==NULL)
      {
         LATL::Print<real_t>(m,n,A,m);
      }
      else
      {
         std::ofstream out(outfile);
         if(out)
            LATL::Print<real_t>(m,n,A,m,out);
         else
         {
            cerr << "Cannot open output file." << endl;
            return 1;
         }
      }
   }
   return 0;
}
