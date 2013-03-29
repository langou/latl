//
//  lagtm.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 6/26/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lagtm_h
#define _lagtm_h

/// @file lagtm.h Performs multiplication by a tridiagonal matrix.


#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Performs multiplication by a real tridiagonal matrix.
   ///
   /// For a tridiagonal matrix A, matrices X and B, and scalars alpha and beta,
   ///
   ///        B := alpha*A*X+beta*B  or  B := alpha*A'*X+beta*B
   /// is computed.  The scalars alpha and beta must be either 1,-1 or 0.  
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies if the transpose of A is used:
   ///
   ///        'N' or 'n': B := alpha*A*X+beta*B
   ///        'T' or 't': B := alpha*A'*X+beta*B
   ///        'C' or 'c': B := alpha*A'*X+beta*B   
   /// @param n The order of the tridiagonal matrix A.  n>=0
   /// @param nrhs The number of columns of matrices X and B.  nrhs>=0
   /// @param alpha Real scalar, must be either 1, -1 or 0; otherwise, assumed to be 0.
   /// @param dl Pointer to vector of length n-1 containing the subdiagonal of A.
   /// @param d Pointer to vector of length n containing the diagonal of A.
   /// @param du Pointer to vector of length n-1 containing the superdiagonal of A.
   /// @param X Pointer to n-by-nrhs matrix X.
   /// @param ldX Column length of X.  ldX>=n
   /// @param beta Real scalar, must be either 1, -1 or 0; otherwise, assumed to be 1.
   /// @param B Pointer to n-by-nrhs matrix B.  On exit, B contains the result of the computation.
   /// @param ldB Column length of B.  ldB>=n
   /// @ingroup AUX
   
   template<typename real_t>
   int LAGTM(char trans,int_t n,int_t nrhs,real_t alpha,real_t *dl,real_t *d,real_t *du,real_t *X,int_t ldX,real_t beta,real_t *B,int_t ldB)
   {
      using std::toupper;
      trans=toupper(trans);
      if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -1;
      else if(n<0)
         return -2;
      else if(nrhs<0)
         return -3;
      else if(ldX<n)
         return -9;
      else if(ldB<n)
         return -12;
      
      const real_t one(1.0);
      const real_t zero(0.0);
      
      if(n>0)
      {
         if(beta==zero)
         {
            real_t *b=B;
            for(int_t j=0;j<nrhs;j++)
            {
               for(int_t i=0;i<n;i++)
                  b[i]=zero;
               b+=ldB;
            }
         }
         else if(beta==-one)
         {
            real_t *b=B;
            for(int_t j=0;j<nrhs;j++)
            {
               for(int_t i=0;i<n;i++)
                  b[i]=-b[i];
               b+=ldB;
            }
         }
         if(alpha==one)
         {
            if(trans=='N')
            {
               real_t *b=B;
               real_t *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=d[0]*x[0];
                  }
                  else
                  {
                     b[0]+=d[0]*x[0]+du[0]*x[1];
                     b[n-1]+=dl[n-2]*x[n-2]+d[n-1]*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=dl[i-1]*x[i-1]+d[i]*x[i]+du[i]*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
            else
            {
               real_t *b=B;
               real_t *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=d[0]*x[0];
                  }
                  else
                  {
                     b[0]+=d[0]*x[0]+dl[0]*x[1];
                     b[n-1]+=du[n-2]*x[n-2]+d[n-1]*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=du[i-1]*x[i-1]+d[i]*x[i]+dl[i]*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
         }
         else if(alpha==-one)
         {
            if(trans=='N')
            {
               real_t *b=B;
               real_t *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=-d[0]*x[0];
                  }
                  else
                  {
                     b[0]+=-d[0]*x[0]-du[0]*x[1];
                     b[n-1]+=-dl[n-2]*x[n-2]-d[n-1]*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=-dl[i-1]*x[i-1]-d[i]*x[i]-du[i]*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
            else
            {
               real_t *b=B;
               real_t *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=-d[0]*x[0];
                  }
                  else
                  {
                     b[0]+=-d[0]*x[0]-dl[0]*x[1];
                     b[n-1]+=-du[n-2]*x[n-2]-d[n-1]*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=-du[i-1]*x[i-1]-d[i]*x[i]-dl[i]*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
         }
      }
      return info;
   }

   /// @brief Performs multiplication by a complex tridiagonal matrix.
   ///
   /// For a tridiagonal matrix A, matrices X and B, and scalars alpha and beta, one of
   ///
   ///        B := alpha*A*X+beta*B
   ///        B := alpha*A.'*X+beta*B
   ///        B := alpha*A'*X+beta*B
   /// is computed.  The scalars alpha and beta must be either 1,-1 or 0.  
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param trans Specifies if the transpose or conjugate transpose of A is used:
   ///
   ///        'N' or 'n': B := alpha*A*X+beta*B
   ///        'T' or 't': B := alpha*A.'*X+beta*B
   ///        'C' or 'c': B := alpha*A'*X+beta*B   
   /// @param n The order of the tridiagonal matrix A.  n>=0
   /// @param nrhs The number of columns of matrices X and B.  nrhs>=0
   /// @param alpha Real scalar, must be either 1, -1 or 0; otherwise, assumed to be 0.
   /// @param dl Pointer to vector of length n-1 containing the subdiagonal of A.
   /// @param d Pointer to vector of length n containing the diagonal of A.
   /// @param du Pointer to vector of length n-1 containing the superdiagonal of A.
   /// @param X Pointer to n-by-nrhs matrix X.
   /// @param ldX Column length of X.  ldX>=n
   /// @param beta Real scalar, must be either 1, -1 or 0; otherwise, assumed to be 1.
   /// @param B Pointer to n-by-nrhs matrix B.  On exit, B contains the result of the computation.
   /// @param ldB Column length of B.  ldB>=n
   /// @ingroup AUX
   
   template<typename real_t>
   int LAGTM(char trans,int_t n,int_t nrhs,real_t alpha,complex<real_t> *dl,complex<real_t> *d,complex<real_t> *du,complex<real_t> *X,int_t ldX,real_t beta,complex<real_t> *B,int_t ldB)
   {
      using std::conj;
      using std::toupper;
      trans=toupper(trans);
      if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -1;
      else if(n<0)
         return -2;
      else if(nrhs<0)
         return -3;
      else if(ldX<n)
         return -9;
      else if(ldB<n)
         return -12;
      
      const real_t one(1.0);
      const real_t zero(0.0);
      
      if(n>0)
      {
         if(beta==zero)
         {
            complex<real_t> *b=B;
            for(int_t j=0;j<nrhs;j++)
            {
               for(int_t i=0;i<n;i++)
                  b[i]=zero;
               b+=ldB;
            }
         }
         else if(beta==-one)
         {
            complex<real_t> *b=B;
            for(int_t j=0;j<nrhs;j++)
            {
               for(int_t i=0;i<n;i++)
                  b[i]=-b[i];
               b+=ldB;
            }
         }
         if(alpha==one)
         {
            if(trans=='N')
            {
               complex<real_t> *b=B;
               complex<real_t> *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=d[0]*x[0];
                  }
                  else
                  {
                     b[0]+=d[0]*x[0]+du[0]*x[1];
                     b[n-1]+=dl[n-2]*x[n-2]+d[n-1]*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=dl[i-1]*x[i-1]+d[i]*x[i]+du[i]*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
            else if(trans=='T')
            {
               complex<real_t> *b=B;
               complex<real_t> *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=d[0]*x[0];
                  }
                  else
                  {
                     b[0]+=d[0]*x[0]+dl[0]*x[1];
                     b[n-1]+=du[n-2]*x[n-2]+d[n-1]*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=du[i-1]*x[i-1]+d[i]*x[i]+dl[i]*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
            else if(trans=='C')
            {
               complex<real_t> *b=B;
               complex<real_t> *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=conj(d[0])*x[0];
                  }
                  else
                  {
                     b[0]+=conj(d[0])*x[0]+conj(dl[0])*x[1];
                     b[n-1]+=conj(du[n-2])*x[n-2]+conj(d[n-1])*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=conj(du[i-1])*x[i-1]+conj(d[i])*x[i]+conj(dl[i])*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
         }
         else if(alpha==-one)
         {
            if(trans=='N')
            {
               complex<real_t> *b=B;
               complex<real_t> *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=-d[0]*x[0];
                  }
                  else
                  {
                     b[0]+=-d[0]*x[0]-du[0]*x[1];
                     b[n-1]+=-dl[n-2]*x[n-2]-d[n-1]*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=-dl[i-1]*x[i-1]-d[i]*x[i]-du[i]*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
            else if(trans=='T')
            {
               complex<real_t> *b=B;
               complex<real_t> *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=-d[0]*x[0];
                  }
                  else
                  {
                     b[0]+=-d[0]*x[0]-dl[0]*x[1];
                     b[n-1]+=-du[n-2]*x[n-2]-d[n-1]*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=-du[i-1]*x[i-1]-d[i]*x[i]-dl[i]*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
            else if(trans=='C')
            {
               complex<real_t> *b=B;
               complex<real_t> *x=X;
               for(int_t j=0;j<nrhs;j++)
               {
                  if(n==1)
                  {
                     b[0]+=-conj(d[0])*x[0];
                  }
                  else
                  {
                     b[0]+=-conj(d[0])*x[0]-conj(dl[0])*x[1];
                     b[n-1]+=-conj(du[n-2])*x[n-2]-conj(d[n-1])*x[n-1];
                     for(int_t i=1;i<n-1;i++)
                        b[i]+=-conj(du[i-1])*x[i-1]-conj(d[i])*x[i]-conj(dl[i])*x[i+1];
                  }
                  b+=ldB;
                  x+=ldX;
               }
            }
         }
      }
      return info;
   }
}
#endif
