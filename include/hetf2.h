//
//  hetf2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 6/26/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _hetf2_h
#define _hetf2_h

/// @file hetf2.h Computes the factorization of a complex Hermitian matrix A using the Bunch-Kaufman diagonal pivoting method.

#include <cmath>
#include "imax.h"
#include "lapy2.h"
#include "swap.h"
#include "her.h"
#include "scal.h"
#include "latl.h"

namespace latl
{
   /// @brief Computes the factorization of a complex Hermitian matrix A using the Bunch-Kaufman diagonal pivoting method:
   ///
   ///        A = U * D * U^H     or  A = L * D * L^H
   /// where U or L is a product of permutation and unit upper or lower triangular matrices.
   /// U^H is the conjugate transpose of U, and D is Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith column has a zero pivot with no other available pivots, or is NaN
   /// @tparam real_t Floating point type.
   /// @param uplo Input character:
   ///
   ///     'U' or 'u' designates that only the upper triangular portion of the matrix A will be accessed
   ///     'L' or 'l' designates that only the lower triangular portion of the matrix A will be accessed
   /// @param n The order of the matrix A.  n >= 0
   /// @param A Complex matrix size ldA-by-n. On entry, the Hermitian matrix A.  On exit, the block diagonal matrix D and the multipliers used to obtain the factor U or L.
   /// @param ldA The column length of the array A.
   /// @param IPIV Integer array size n.  On exit, contains the details of the interchanges and the block structure of D.
   /// @param BSDV Boolean array size n.
   /// If uplo=='U' and BSDV[k] and BSDV[k-1] are true, then rows and
   /// columns k-1 and IPIV[k] were interchanged and D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
   /// If uplo=='L' and BSDV[k] and BSDV[k+1] are true, then rows and columns k+1 and IPIV[k]
   /// were interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.

   template <typename real_t>
   int_t hetf2(const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA, int_t * const IPIV, bool * const BSDV)
   {
      if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      const real_t alpha = (1.0 + std::sqrt(17.0))/8.0;
      const real_t one(1.0);
      const real_t zero(0.0);
      int_t info = 0;
      complex<real_t> * Ak;
      int_t k, kstep, imax = 0, kp;
      real_t colmax = 0, absakk, r1;
      
      for (int_t i = 0; i < n; ++i)
         BSDV[i] = 0;
      
      if (uplo == 'U' || uplo == 'u')
      {
         k = n-1,
         kstep = 1;
         Ak = A+ldA*k;
         while (k >= 0)
         {
            kstep = 1;
            Ak = A+ldA*k;
            
            absakk = std::abs(real(Ak[k]));
            
            if (k > 0)
            {
               imax = latl::imax(k, Ak, 1);
               colmax = std::abs(Ak[imax]);
            }
            else
            {
               colmax = zero;
            }
            
            if (std::max(absakk, colmax) == 0 || std::isnan(absakk))
            {
               if (info == 0)
               {
                  info = k+1;
               }
               kp = k;
               Ak[k] = real(Ak[k]);
            }
            else
            {
               if (absakk >= (alpha*colmax))
               {
                  kp = k;
               }
               else
               {
                  int_t jmax;
                  real_t rowmax;
                  
                  jmax = imax + latl::imax(k-imax, A+ldA*(imax+1)+imax, ldA);
                  rowmax = std::abs(*(A+ldA*jmax+imax));
                  if (imax > 0)
                  {
                     jmax = latl::imax(imax-1, A+ldA*imax, 1);
                     rowmax = std::max(rowmax, std::abs(*(A+ldA*imax+jmax)));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                  {
                     kp = k;
                  }
                  else if (std::abs(real(*(A+ldA*imax+imax))) >= (alpha*rowmax))
                  {
                     kp = imax;
                  }
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               
               int_t kinc = k - kstep +1;
               complex<real_t> * Akinc = A + ldA*kinc;
               complex<real_t> * Aj;
               complex<real_t> temp;
               if (kp != kinc)
               {
                  latl::swap(kp-1, Akinc, 1, A+ldA*kp, 1);
                  for (int_t j = kp+1; j< kinc-1; ++j)
                  {
                     Aj = A+ldA*j;
                     temp = conj(Akinc[j]);
                     Akinc[j] = conj(Aj[kp]);
                     Aj[kp] = temp;
                  }
                  Aj = A+ldA*kp;
                  Akinc[kp] = conj(Akinc[kp]);
                  r1 = real(Akinc[kinc]);
                  Akinc[kinc] = real(Aj[kp]);
                  Aj[kp] = r1;
                  if (kstep == 2)
                  {
                     Ak[k] = real(Ak[k]);
                     temp = Ak[k-1];
                     Ak[k-1] = Ak[kp];
                     Ak[kp] = temp;
                  }
               }
               else
               {
                  Ak[k] = real(Ak[k]);
                  if (kstep == 2)
                  {
                     Aj = A+ldA*(k-1);
                     Aj[k-1] = real(Aj[k-1]);
                  }
               }
               
               if (kstep == 1)
               {
                  r1 = one/real(Ak[k]);
                  latl::her(uplo, k, -r1, Ak, 1, A, ldA);
                  latl::scal(k, r1, Ak, 1);
               }
               else
               {
                  if (k > 1)
                  {
                     complex<real_t> * Akm1 = A+ldA*(k-1);
                     // this is the norm of the element above k
                     real_t d = latl::lapy2(real(Ak[k-1]),imag(Ak[k-1]));
                     real_t d22 = real(Akm1[k-1])/d;
                     real_t d11 = real(Ak[k])/d;
                     real_t tt = one/(d11*d22-one);
                     complex<real_t> d12 = Ak[k-1]/d;
                     complex<real_t> wk, wkm1;
                     d = tt/d;
                     for (int_t j = k-2; j >= 0; --j)
                     {
                        Aj = A + ldA*j;
                        wkm1 = d*(d11*Akm1[j]-conj(d12)*Ak[j]);
                        wk = d*(d22*Ak[j]-d12*Akm1[j]);
                        for (int_t i = j; i >= 0; --i)
                        {
                           Aj[i] = Aj[i] - Ak[i]*conj(wk)- Akm1[i]*conj(wkm1);
                        }
                        Ak[j] = wk;
                        Akm1[j] = wkm1;
                        Aj[j] = complex<real_t>(real(Aj[j]), 0.0);
                     }
                  }
               }
            }
            if (kstep == 1)
            {
               IPIV[k] = kp;
            }
            else
            {
               IPIV[k] = kp;
               IPIV[k-1] = kp;
               BSDV[k] = 1;
               BSDV[k-1] = 1;
            }
            
            k = k-kstep;
         } // ends at 435
      }
      else
      {
         k = 0;
         kstep = 1;
         while (k < n)
         {
            kstep = 1;
            Ak = A+ldA*k;
            
            absakk = std::abs(Ak[k]);
            
            if (k < n-1)
            {
               imax = k+1 + latl::imax(n-k-1, Ak+k+1, 1);
               colmax = std::abs(Ak[imax]);
            }
            else
            {
               colmax = zero;
            }
            
            if (std::max(absakk, colmax) == 0 || std::isnan(absakk))
            {
               if (info == 0)
                  info = k+1;
               kp = k;
               Ak[k] = real(Ak[k]);
            }
            else
            {
               if (absakk >= alpha*colmax)
                  kp = k;
               else
               {
                  int_t jmax;
                  real_t rowmax;
                  
                  jmax = k+latl::imax(imax-k, Ak+imax, ldA);
                  rowmax = std::abs(*(A+ldA*jmax+imax));
                  complex<real_t> * Aimax = A+ldA*imax;
                  if (imax < n-1)
                  {
                     jmax = imax + latl::imax(n-imax-1, Aimax+imax+1, 1);
                     rowmax = std::max(rowmax, std::abs(Aimax[jmax]));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                     kp = k;
                  else if (std::abs(real(Aimax[imax])) >= alpha*rowmax)
                  {
                     kp = imax;
                  }
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               
               int_t kinc = k + kstep - 1;
               complex<real_t> * Akinc = A + ldA*kinc;
               complex<real_t> * Aj, * Akp1 = A+ldA*(k+1);
               complex<real_t> temp;
               
               if (kp != kinc)
               {
                  if (kp < n-1)
                  {
                     latl::swap(n-kp-1, Akinc+kp+1, 1, A+ldA*kp+kp+1, 1);
                  }
                  for (int_t j = kinc+1; j < kp; ++j)
                  {
                     Aj = A+ldA*j;
                     temp = conj(Akinc[j]);
                     Akinc[j] = conj(Aj[kp]);
                     Aj[kp] = temp;
                  }
                  Aj = A+ldA*kp;
                  Akinc[kp] = conj(Akinc[kp]);
                  r1 = real(Akinc[kinc]);
                  Akinc[kinc] = real(Aj[kp]);
                  Aj[kp] = r1;
                  if (kstep == 2)
                  {
                     Ak[k] = real(Ak[k]);
                     temp = Ak[k+1];
                     Ak[k+1] = Ak[kp];
                     Ak[kp] = temp;
                  }
               }
               else
               {
                  Ak[k] = real(Ak[k]);
                  if (kstep == 2)
                     Akp1[k+1] = real(Akp1[k+1]);
               }
               if (kstep == 1)
               {
                  if (k < n-1)
                  {
                     r1 = one/real(Ak[k]);
                     latl::her(uplo, n-k-1, -r1, Ak+k+1, 1, Akp1+k+1, ldA);
                     latl::scal(n-k-1, r1, Ak+k+1, 1);
                  }
               }
               else
               {
                  if (k < n-2)
                  {
                     real_t d = latl::lapy2(real(Ak[k+1]),imag(Ak[k+1]));
                     real_t d11 = real(Akp1[k+1])/d;
                     real_t d22 = real(Ak[k])/d;
                     real_t tt = one/(d11*d22-one);
                     complex<real_t> d21 = Ak[k+1]/d;
                     complex<real_t> wk, wkp1;
                     d = tt/d;
                     
                     for (int_t j = k+2; j < n; ++j)
                     {
                        Aj = A + ldA*j;
                        wk = d*(d11*Ak[j]-d21*Akp1[j]);
                        wkp1 = d*(d22*Akp1[j]-conj(d21)*Ak[j]);
                        for (int_t i = j; i < n; ++i)
                        {
                           Aj[i] = Aj[i] - Ak[i]*conj(wk)- Akp1[i]*conj(wkp1);
                        }
                        Ak[j] = wk;
                        Akp1[j] = wkp1;
                        Aj[j] = complex<real_t>(real(Aj[j]), 0.0);
                     }
                  }
               }
               if (kstep == 1)
               {
                  IPIV[k] = kp;
               }
               else
               {
                  IPIV[k] = kp;
                  IPIV[k+1] = kp;
                  BSDV[k] = 1;
                  BSDV[k+1] = 1;
               }
            }
            k += kstep;
         }
      }
      return info;
   }
}


#endif
