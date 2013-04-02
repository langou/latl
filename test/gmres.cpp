#include <iostream>
#include <cstdlib>
#include <nrm2.h>
#include <laset.h>
#include <gemv.h>
#include <copy.h>
#include <rotg.h>
#include <dot.h>
#include <getrf.h>
#include <getri.h>

//
// set floating point type to define Real typedef
//

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

using namespace LATL;
using namespace std;

Real GMRES(const int n, const int m, Real *A, const int ldA, Real *x, Real *b, const Real tol)
{
   const Real zero(0.0);
   const Real one(1.0);
   const int max_iter=n;
   Real *V=new Real[n*(m+1)];
   Real *H=new Real[m*(m+1)];
   Real *cs=new Real[m];
   Real *sn=new Real[m];
   Real *e1=new Real[n];
   Real *s=new Real[n];
   Real *w=new Real[n];
   Real *y=new Real[m];
   e1[0]=one;

   for(int i=0;i<n;i++)
      b[i]=x[i];


   Real bnorm=NRM2(n,b,1);
   if(bnorm==zero)
      bnorm=one;

   Real *r=new Real[n];

   // r := b - Ax

   COPY(n,b,1,r,1);
   GEMV('N',n,n,-one,A,n,x,1,one,r,1);
   Real error=NRM2(n,r,1);
   int iter=0;

   while((iter<max_iter)&&(error>tol))
   {
      COPY(n,b,1,r,1);
      GEMV('N',n,n,-one,A,n,x,1,one,r,1);
      Real temp=NRM2(n,r,1);
      for(int i=0;i<n;i++)
         V[i]=r[i]/temp;
      for(int i=0;i<n;i++)
         s[i]=temp*e1[i];
      Real *v=V;                                       // v = V(:,i)
      Real *h=H;                                       // h = H(:,i)
      for(int i=0;i<m;i++)
      {
         GEMV('N',n,n,one,A,n,v,1,zero,w,1);           // w := Av
         Real *u=V;
         for(int k=0;k<=i;k++)
         {
            h[k]=DOT(n,w,1,u,1);
            for(int j=0;j<n;j++)
               w[j]-=h[k]*u[j];
            u+=n;
         }
         h[i+1]=NRM2(n,w,1);
         for(int j=0;j<n;j++)
            u[j]=w[j]/h[i+1];
         for(int k=0;k<i;k++)
         {
            temp =  cs[k]*h[k] + sn[k]*h[k+1];
            h[k+1]    = -sn[k]*h[k] + cs[k]*h[k+1];
            h[k]      =  temp;
         }
         ROTG(h[i],h[i+1],cs[i],sn[i]);
         temp=cs[i]*s[i];
         s[i+1]=-sn[i]*s[i];
         s[i]=temp;
         error=abs(s[i+1])/bnorm;
         if(error>tol)
         {
            int l=i+1;
            int_t *ipiv=new int_t[l];
            Real *M=new Real[l*l];
            for(int J=0;J<l;J++)
               for(int I=0;I<l;I++)
                  M[I+J*l]=H[I+J*(m+1)];
            GETRF(l,l,M,l,ipiv);
            GETRI(l,M,l,ipiv,l);
            GEMV('N',l,l,one,M,l,s,1,zero,y,1);
            GEMV('N',n,l,one,V,n,y,1,one,x,1);
         }
         v+=n; // update v to next column of V
         h+=m+1; // update h to next column of H
      }
      if(error>tol)
      {
         int_t *ipiv=new int_t[m];
         Real *M=new Real[m*m];
         for(int J=0;J<m;J++)
            for(int I=0;I<m;I++)
               M[I+J*m]=H[I+J*(m+1)];
         GETRF(m,m,M,m,ipiv);
         GETRI(m,M,m,ipiv,m);
         GEMV('N',m,m,one,M,m,s,1,zero,y,1);
         GEMV('N',n,m,one,V,n,y,1,one,x,1);
         COPY(n,b,1,r,1);
         GEMV('N',n,n,-one,A,n,x,1,one,r,1);
         error=NRM2(n,r,1)/bnorm;
      }
   }

   delete [] y;
   delete [] w;
   delete [] b;
   delete [] e1;
   delete [] sn;
   delete [] cs;
   delete [] H;
   delete [] V;
   
   return error;
}

int main(int argc, char **argv)
{
   int n=100;
   if(argc>1)
      n=atoi(argv[1]);
   Real *A=new Real[n*n];
   Real *b=new Real[n];
   Real *x=new Real[n];
   Real tol=1.0e-07;

   for(int i=0;i<n;i++)
   {
      x[i]=1.0;
      b[i]=static_cast<Real>(i)/static_cast<Real>(n);
      for(int j=0;j<n;j++)
      {
         if(i==j)
            A[i+j*n]=2.0;
         else if(j==i-1)
            A[i+j*n]=1.0;
         else if(j==i+1)
            A[i+j*n]=-1.0;
      }
   }

   Real error=GMRES(n,n/2,A,n,x,b,tol);

   cerr << "error=" << error << endl;
   for(int i=0;i<n;i++)
      cout << x[i] << endl;
   
   return 0;
}