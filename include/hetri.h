//
//  hetri.h
//  Linear Algebra Template Library
//
//  Created by Henricus Bouwmeester on 1/25/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _hetri_h
#define _hetri_h

/// @file hetri.h Computes the inverse of a Hermitian indefinite matrix.

#include "latl.h"
#include "dot.h"
#include "copy.h"
#include "swap.h"
#include "hemv.h"

namespace latl
{
   /// @brief Computes the inverse of a Hermitian indefinite matrix.
   /// 

   template <typename real_t>
      int_t hetri(char uplo, int_t n, complex<real_t> *A, int_t ldA, int_t *ipiv, int_t nb=32)
      {

         using std::toupper;
         uplo=toupper(uplo);
         const complex<real_t> zero = complex<double>(0.0,0.0);
         const complex<real_t> one = complex<double>(1.0,0.0);

         if((uplo!='U')&&(uplo!='L'))
            return -1;
         else if(n<0)
            return -2;
         else if(ldA<n)
            return -4;
         else if(n==0)
            return 0;

         if(uplo=='U')
         {
            for(int_t i=0; i<n; i++)
               if((ipiv[i] > 0) && (A[i+i*ldA]==zero))
                  return 1;
         }
         else
         {
            for(int_t i=n-1; i>=0; i--)
               if((ipiv[i] > 0) && (A[i+i*ldA]==zero))
                  return 1;
         }


         complex<real_t> *work =new complex<real_t>[n];
         int_t kstep;

         if (uplo=='U')
         {
            int_t k = 1;
            while(k<=n)
            {
               if(ipiv[k-1]>-1)
               {
                  A[(k-1)+(k-1)*ldA].real() = real(one) / real(A[(k-1)+(k-1)*ldA]);

                  if(k>1)
                  {
                     copy<real_t>( k-1, A+(k-1)*ldA, 1, work, 1);
                     hemv<real_t>( uplo, k-1, -one, A, ldA, work, 1, zero, A+(k-1)*ldA, 1);
                     A[(k-1)+(k-1)*ldA].real() -= real(dotc( k-1, work, 1, A+(k-1)*ldA, 1));
                  }
                  kstep = 1;
               }
               else
               {
                  complex<real_t> t;
                  complex<real_t> Ak;
                  complex<real_t> Akp1;
                  complex<real_t> Akkp1;
                  complex<real_t> d;
                  t = std::abs( A[(k-1)+k*ldA] );
                  Ak.real() = real(A[(k-1)+(k-1)*ldA]) / real(t);
                  Akp1.real() = real(A[k+k*ldA]) / real(t);
                  Akkp1 = A[(k-1)+k*ldA] / t;
                  d.real() = real(t)*( real(Ak)*real(Akp1) - real(one) );
                  std::cout << "d = " << d << std::endl;
                  A[(k-1)+(k-1)*ldA].real() = real(Akp1) / real(d);
                  A[k+k*ldA].real() = real(Ak) / real(d);
                  A[(k-1)+k*ldA] = Akkp1 / d;

                  if(k>1)
                  {
                     copy<real_t>( k-1, A+(k-1)*ldA, 1, work, 1);
                     hemv<real_t>( uplo, k-1, -one, A, ldA, work, 1, zero, A+(k-1)*ldA, 1);
                     A[(k-1)+(k-1)*ldA].real() -= real(dotc( k-1, work, 1, A+(k-1)*ldA, 1));
                     A[(k-1)+k*ldA] -= dotc( k-1, A+(k-1)*ldA, 1, A+k*ldA, 1);
                     copy<real_t>( k-1, A+k*ldA, 1, work, 1);
                     hemv<real_t>( uplo, k-1, -one, A, ldA, work, 1, zero, A+k*ldA, 1);
                     A[k+k*ldA].real() -= real(dotc( k-1, work, 1, A+k*ldA, 1));
                  }
                  kstep = 2;
               }
               
               int_t kp = std::abs( ipiv[k-1] );
//                for (int i=0; i<n; i++)
//                std::cout << "ipiv[" << i << "] = " << ipiv[i] << std::endl;
return 0;
               if(kp != k-1)
               {
                  std::cout << "k = " << k << ", kp = " << kp << std::endl;
                  swap<real_t>( kp, A+(k-1)*ldA, 1, A+kp*ldA, 1);
return 0;
                  for(int_t j=kp+2;j<=k-1;j++)
                  {
                     complex<real_t> temp = conj(A[(j-1)+(k-1)*ldA]);
                     A[(j-1)+(k-1)*ldA] = conj(A[kp+(j-1)*ldA]);
                     A[kp+(j-1)*ldA] = temp;
                  }

                  A[kp+(k-1)*ldA] = conj(A[kp+(k-1)*ldA]);
                  complex<real_t> temp = A[(k-1)+(k-1)*ldA];
                  A[(k-1)+(k-1)*ldA] = A[kp+kp*ldA];
                  A[kp+kp*ldA] = temp;
                  if(kstep==2)
                  {
                     temp = A[(k-1)+k*ldA];
                     A[(k-1)+k*ldA] = A[kp+k*ldA];
                     A[kp+k*ldA] = temp;
                  }
               }
               k = k + kstep;
            }
         }

         return 0;
      }
}

#endif
