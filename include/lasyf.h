//
//  lasyf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 7/25/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _lasyf_h
#define _lasyf_h

/// @file lasyf.h Computes a partial factorization of a real symmetric matrix A.

#include "copy.h"
#include "gemm.h"
#include "gemv.h"
#include "imax.h"
#include "scal.h"
#include "swap.h"
#include "latl.h"

namespace LATL
{
   /// @brief Computes a partial factorization of a real symmetric matrix A using the Bunch-Kaufman diagonal pivoting method.  The partial factorization has the form:
   ///
   ///     A = (I  U12 ) (A11 0 ) ( I   0   )   if uplo = 'U'
   ///         (0  U22 ) (0   D ) (U12' U22')
   ///     A = (L11  0 ) (D   0 ) (L11' L21')   if uplo = 'L'
   ///         (L21  I ) (0  A22) ( 0    I  )
   ///
   /// where the order of D is at most nb.  The actual order is returned in the argument kb, and is either nb or nb-1, or n if n <= nb.
   /// This routine is an auxiliary routine for sytrf.  It uses blocked code to update the submatrix A11 (if uplo = 'U') or A22 (if uplo = 'L').
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith pivot of D is exactly zero.  The factorization has been completed, but the block diagonal matrix D is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetric matrix A is stored.
   /// @param n Order of the matrix A.  n >= 0.
   /// @param nb The maximum number of columns of the matrix A that should be factored.  nb should be at least 2 to allow for 2-by-2 pivot blocks.
   /// @param kb The number of columns of A that were actually factored.  kb is either nb-1 or nb, or n if n <= nb.
   /// @param A Real array size ldA-by-n.  On entry, the symmetric matrix A.   On exit, A contains details of the partial factorization.
   /// @param ldA Column length of the array A.  ldA >= n
   /// @param ipiv Integer array size n. On exit, details of the interchanges of D. If uplo = 'U', the last kb elements of ipiv and bsdv are set, and if uplo = 'L', the first kb elements of ipiv and bsdv are set.
   /// @param bsdv Bool array size n. On exit, contains the details of the block structure of D.  If bsdv[k] = 0, then rows and columns k and ipiv[k] were interchanged and D[k, k] is a 1-by-1 diagonal block.  If bsdv[k] = 1, then k is part of a 2-by-2 diagonal block.  In a 2 by 2 block, if uplo = 'U', and ipiv[k] = ipiv[k-1], then rows and columns k-1 and ipiv[k] were interchanged.  If uplo = 'L' and ipiv[k] = ipiv[k+1], then rows and columns k+1 and ipiv[k] were interchanged.
   /// @param Work Real array size n-by-nb.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t LASYF(const char uplo, const int_t n, const int_t nb, int_t &kb, real_t * const A, const int_t ldA, int_t * const ipiv, bool * const bsdv, real_t * Work)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (nb < 2)
         return -3;
      if (ldA < n)
         return -6;
      
      if (n == 0)
         return 0;
      
      using std::abs;
      using std::min;
      using std::max;
      const real_t alpha = (1.0+sqrt(17.0))/8.0;
      const real_t one(1.0);
      const real_t zero(0.0);
      const int_t ldWork = n;
      int_t k, kw, kp, kstep, imax, jmax, kk, kkw, jtemp, jb, jp;
      real_t * Ak, absakk, colmax, *Aimax, rowmax, *Akk, r1, d21, d11, d22, *Akm1;
      int_t info = 0;
      
      if (uplo == 'U' || uplo == 'u')
      {
         k = n-1;
         int_t kp1;
         real_t *Wkwm1, *Wkw, *Wkkw;
         while (k >= 0)
         {
            kw = nb+k-n;
            Ak = A+ldA*k;
            Akm1 = Ak-ldA;
            Wkw = Work+ldWork*kw;
            Wkwm1 = Wkw-ldWork;
            
            if (nb < n && k <= n-nb )
               break;
            LATL::COPY(k+1, Ak, 1, Wkw, 1);
            if (k < n-1)
            {
               LATL::GEMV('N', k+1, n-k-1, -one, Ak+ldA, ldA, Wkw+ldWork+k, ldWork, one, Wkw, 1);
            }
            
            kstep = 1;
            
            absakk = abs(Wkw[k]);
            
            if (k > 0)
            {
               imax = LATL::IMAX(k, Wkw, 1);
               colmax = abs(Wkw[imax]);
            }
            else
            {
               colmax = zero;
            }
            
            if (max(absakk, colmax) == zero)
            {
               if (info == 0)
               {
                  info = k+1;
               }
               kp = k;
            }
            else
            {
               if (absakk >= alpha*colmax)
               {
                  kp = k;
               }
               else
               {
                  Aimax = A + ldA*imax;
                  LATL::COPY(imax +1, Aimax, 1, Wkwm1, 1);
                  LATL::COPY(k-imax, Aimax+ldA+imax, ldA, Wkwm1+imax+1, 1);
                  
                  if (k < n-1)
                  {
                     LATL::GEMV('N', k+1, n-k-1, -one, Ak+ldA, ldA, Wkw+ldWork+imax, ldWork, one, Wkwm1, 1);
                  }
                  
                  jmax = imax+1 + LATL::IMAX(k-imax, Wkwm1+imax+1, 1);
                  rowmax = abs(Wkwm1[jmax]);
                  if (imax > 0)
                  {
                     jmax = LATL::IMAX(imax, Wkwm1, 1);
                     rowmax = max(rowmax, abs(Wkwm1[jmax]));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                     kp = k;
                  else if (abs(Wkwm1[imax]) >= alpha*rowmax)
                  {
                     kp = imax;
                     LATL::COPY(k+1, Wkwm1, 1, Wkw, 1);
                  }
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               
               kk = k - kstep +1;
               kkw = nb + kk -n;
               Akk = A + ldA*kk;
               Wkkw = Work + ldWork*kkw;
               if (kp != kk)
               {
                  Ak[kp] = Ak[kk];
                  LATL::COPY(k-1-kp, Akk+kp+1, 1, A + ldA*(kp+1)+ kp, ldA);
                  LATL::COPY(kp+1, Akk, 1, A+ldA*kp, 1);
                  LATL::SWAP(n-kk, Akk+kk, ldA, Akk+kp, ldA);
                  LATL::SWAP(n-kk, Wkkw+kk, ldWork, Wkkw+kp, ldWork);
               }
               
               if (kstep == 1)
               {
                  LATL::COPY(k+1, Wkw, 1, Ak, 1);
                  r1 = one/Ak[k];
                  LATL::SCAL(k, r1, Ak, 1);
               }
               else
               {
                  if (k > 1)
                  {
                     d21 = Wkw[k-1];
                     d11 = Wkw[k]/d21;
                     d22 = Wkwm1[k-1]/d21;
                     d21 = (one/(d11*d22-one))/d21;
                     for (int_t j = 0; j <= k-2; ++j)
                     {
                        Akm1[j] = d21*(d11*Wkwm1[j]-Wkw[j]);
                        Ak[j] = d21*(d22*Wkw[j] - Wkwm1[j]);
                     }
                  }
                  
                  Akm1[k-1] = Wkwm1[k-1];
                  Ak[k-1] = Wkw[k-1];
                  Ak[k] = Wkw[k];
               }
               
            }
            if (kstep == 1)
            {
               ipiv[k] = kp;
               bsdv[k] = 0;
            }
            else
            {
               ipiv[k] = kp;
               ipiv[k-1] = kp;
               bsdv[k] = 1;
               bsdv[k-1] = 1;
            }
            
            k -= kstep;
         }
         
         for (int_t j = ((k)/nb)*nb; j >= 0; j -= nb)
         {
            jb = min(nb, k-j+1);
            for (int_t jj = j; jj < j+jb; ++jj)
            {
               LATL::GEMV('N', jj-j+1, n-k-1,-one, Ak+ldA+j, ldA, Wkw+ldWork+jj, ldWork, one, A+jj*ldA+j, 1);
            }
            LATL::GEMM('N', 'T', j, jb, n-k-1, -one, Ak+ldA, ldA, Wkw+ldWork+j, ldWork, one, A+j*ldA, ldA);
         }
         
         kp1 = k+1;
         while (kp1 < n)
         {
            jtemp = kp1;
            jp = ipiv[kp1];
            if (bsdv[kp1] == 1)
            {
               ++kp1;
            }
            ++kp1;
            if (jp != jtemp && kp1 < n)
            {
               LATL::SWAP(n-kp1, A+ldA*kp1+jp, ldA, A+ldA*kp1+jtemp, ldA);
            }
         }
         kb = n-k-1;
         
      }
      else
      {
         k = 0;
         real_t * Wk, *Wkp1, *Akp1;
         int_t km1;
         while (k < n)
         {
            if (k >= nb-1 && nb < n)
               break;
            
            Ak = A+ldA*k;
            Akp1 = Ak+ldA;
            Wk = Work+ldWork*k;
            Wkp1 = Wk+ldWork;
            LATL::COPY(n-k, Ak+k, 1, Wk+k, 1);
            LATL::GEMV('N', n-k, k, -one, A+k, ldA, Work+k, ldWork, one, Wk+k, 1);
            
            kstep = 1;
            
            absakk = abs(Wk[k]);
            
            if (k < n-1)
            {
               imax = k+1 + LATL::IMAX(n-k-1, Wk+k+1, 1);
               colmax = abs(Wk[imax]);
            }
            else
            {
               colmax = zero;
            }
            
            if (max(absakk, colmax) == zero)
            {
               if (info == 0)
                  info = k+1;
               kp = k;
            }
            else
            {
               if (absakk >= alpha*colmax)
                  kp = k;
               else
               {
                  Aimax = A+ldA*imax;
                  LATL::COPY(imax-k, Ak+imax, ldA, Wkp1+k, 1);
                  LATL::COPY(n-imax, Aimax+imax, 1, Wkp1+imax, 1);
                  LATL::GEMV('N', n-k, k, -one, A+k, ldA, Work+imax, ldWork, one, Wkp1+k, 1);
                  
                  jmax = k + LATL::IMAX(imax-k, Wkp1+k, 1);
                  
                  rowmax = abs(Wkp1[jmax]);
                  if (imax < n-1)
                  {
                     jmax = imax+1+LATL::IMAX(n-imax-1, Wkp1+imax+1, 1);
                     rowmax = max(rowmax, abs(Wkp1[jmax]));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                  {
                     kp = k;
                  }
                  else if ( abs(Wkp1[imax]) >= alpha*rowmax)
                  {
                     kp = imax;
                     LATL::COPY(n-k, Wkp1+k, 1, Wk+k, 1);
                  }
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               kk = k+kstep-1;
               Akk = A+ldA*kk;
               
               if (kp != kk)
               {
                  Ak[kp] = Ak[kk];
                  LATL::COPY(kp-k-1, Akk+k+1, 1, Akp1+kp, ldA);
                  LATL::COPY(n-kp, Akk+kp, 1, A+ldA*kp+kp, 1);
                  LATL::SWAP(kk+1, A+kk, ldA, A+kp, ldA);
                  LATL::SWAP(kk+1, Work+kk, ldWork, Work+kp, ldWork);
               }
               
               if (kstep == 1)
               {
                  LATL::COPY(n-k, Wk+k, 1, Ak+k, 1);
                  if (k < n-1)
                  {
                     r1 = one/Ak[k];
                     LATL::SCAL(n-k-1, r1, Ak+k+1, 1);
                  }
               }
               else
               {
                  if ( k < n-2)
                  {
                     d21 = Wk[k+1];
                     d11 = Wkp1[k+1]/d21;
                     d22 = Wk[k]/d21;
                     d21 = (one/(d11*d22-one))/d21;
                     for (int_t j = k+2; j < n; ++j)
                     {
                        Ak[j] = d21*(d11*Wk[j] - Wkp1[j]);
                        Akp1[j] = d21*(d22*Wkp1[j]-Wk[j]);
                     }
                  }
                  
                  Ak[k] = Wk[k];
                  Ak[k+1] = Wk[k+1];
                  Akp1[k+1] = Wkp1[k+1];
               }
            }
            
            if (kstep == 1)
            {
               ipiv[k] = kp;
               bsdv[k] = 0;
            }
            else
            {
               ipiv[k] = kp;
               ipiv[k+1] = kp;
               bsdv[k] = 1;
               bsdv[k+1] = 1;
            }
            
            k += kstep;
         }
         for (int_t j = k; j < n; j += nb)
         {
            jb = min(nb, n-j);
            for (int_t jj = j; jj < j+jb; ++jj)
            {
               LATL::GEMV('N', j+jb-jj, k, -one, A+jj, ldA, Work+jj, ldWork, one, A+jj*ldA+jj, 1);
            }
            
            if (j+jb < n)
            {
               LATL::GEMM('N', 'T', n-j-jb, jb, k, -one, A+j+jb, ldA, Work+j, ldWork, one, A+j*ldA+j+jb, ldA);
            }
         }
         
         km1 = k-1;
         while (km1 >= 0)
         {
            jtemp = km1;
            jp = ipiv[km1];
            if (bsdv[km1] == 1)
            {
               --km1;
            }
            --km1;
            if (jp != jtemp && km1 >= 0)
            {
               LATL::SWAP(km1+1, A+jp, ldA, A+jtemp, ldA);
            }
         }
         
         kb = k;
      }
      
      return info;
   }
   
   /// @brief Computes a partial factorization of a complex symmetric matrix A using the Bunch-Kaufman diagonal pivoting method.  The partial factorization has the form:
   ///
   ///     A = (I  U12 ) (A11 0 ) ( I   0   )   if uplo = 'U'
   ///         (0  U22 ) (0   D ) (U12' U22')
   ///     A = (L11  0 ) (D   0 ) (L11' L21')   if uplo = 'L'
   ///         (L21  I ) (0  A22) ( 0    I  )
   ///
   /// where the order of D is at most nb.  The actual order is returned in the argument kb, and is either nb or nb-1, or n if n <= nb.
   /// This routine is an auxiliary routine for sytrf.  It uses blocked code to update the submatrix A11 (if uplo = 'U') or A22 (if uplo = 'L').
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return i+1 if the ith pivot of D is exactly zero.  The factorization has been completed, but the block diagonal matrix D is exactly singular.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetric matrix A is stored.
   /// @param n Order of the matrix A.  n >= 0.
   /// @param nb The maximum number of columns of the matrix A that should be factored.  nb should be at least 2 to allow for 2-by-2 pivot blocks.
   /// @param kb The number of columns of A that were actually factored.  kb is either nb-1 or nb, or n if n <= nb.
   /// @param A Complex array size ldA-by-n.  On entry, the symmetric matrix A.   On exit, A contains details of the partial factorization.
   /// @param ldA Column length of the array A.  ldA >= n
   /// @param ipiv Integer array size n. On exit, details of the interchanges of D. If uplo = 'U', the last kb elements of ipiv and bsdv are set, and if uplo = 'L', the first kb elements of ipiv and bsdv are set.
   /// @param bsdv Bool array size n. On exit, contains the details of the block structure of D.  If bsdv[k] = 0, then rows and columns k and ipiv[k] were interchanged and D[k, k] is a 1-by-1 diagonal block.  If bsdv[k] = 1, then k is part of a 2-by-2 diagonal block.  In a 2 by 2 block, if uplo = 'U', and ipiv[k] = ipiv[k-1], then rows and columns k-1 and ipiv[k] were interchanged.  If uplo = 'L' and ipiv[k] = ipiv[k+1], then rows and columns k+1 and ipiv[k] were interchanged.
   /// @param Work Complex array size n-by-nb.
   /// @ingroup TRF
   
   template< typename real_t>
   int_t LASYF(const char uplo, const int_t n, const int_t nb, int_t &kb, complex<real_t> * const A, const int_t ldA, int_t * const ipiv, bool * const bsdv, complex<real_t> *Work)
   {
      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (nb < 2)
         return -3;
      if (ldA < n)
         return -6;
      
      if (n == 0)
         return 0;
      
      using std::abs;
      using std::max;
      using std::min;
      const real_t alpha = (1.0+sqrt(17.0))/8.0;
      const real_t zero(0.0);
      const complex<real_t> onec(1.0);
      const int_t ldWork = n;
      int_t k, kw, kp, kstep, imax, jmax, kk, kkw, jtemp, jb, jp;
      real_t absakk, colmax, rowmax;
      complex<real_t> * Ak, *Aimax,  *Akk, *Wkkw, d12, d21, d11, d22, r1, *Akm1;
      int_t info = 0;
      
      if (uplo == 'U' || uplo == 'u')
      {
         k = n-1;
         int_t kp1;
         complex<real_t> *Wkwm1, *Wkw;
         while (k >= 0)
         {
            kw = nb+k-n;
            
            Ak = A+ldA*k;
            Akm1 = Ak-ldA;
            Wkw = Work+ldWork*kw;
            Wkwm1 = Wkw-ldWork;
            
            
            if (nb < n && k <= n-nb )
               break;
            LATL::COPY(k+1, Ak, 1, Wkw, 1);
            if (k < n-1)
            {
               LATL::GEMV('N', k+1, n-k-1, -onec, Ak+ldA, ldA, Wkw+ldWork+k, ldWork, onec, Wkw, 1);
            }
            
            kstep = 1;
            
            absakk = abs(real(Wkw[k]))+abs(imag(Wkw[k]));
            
            if (k > 0)
            {
               imax = LATL::IMAX(k, Wkw, 1);
               colmax = abs(real(Wkw[imax]))+abs(imag(Wkw[imax]));
            }
            else
            {
               colmax = zero;
            }
            
            if (max(absakk, colmax) == zero)
            {
               if (info == 0)
               {
                  info = k+1;
               }
               kp = k;
            }
            else
            {
               if (absakk >= alpha*colmax)
               {
                  kp = k;
               }
               else
               {
                  Aimax = A + ldA*imax;
                  LATL::COPY(imax+1, Aimax, 1, Wkwm1, 1);
                  LATL::COPY(k-imax, Aimax+ldA+imax, ldA, Wkwm1+imax+1, 1);
                  
                  if (k < n-1)
                  {
                     LATL::GEMV('N', k+1, n-k-1, -onec, Ak+ldA, ldA, Wkw+ldWork+imax, ldWork, onec, Wkwm1, 1);
                  }
                  
                  jmax = imax+1 + LATL::IMAX(k-imax, Wkwm1+imax+1, 1);
                  rowmax = abs(real(Wkwm1[jmax]))+abs(imag(Wkwm1[jmax]));
                  if (imax > 0)
                  {
                     jmax = LATL::IMAX(imax, Wkwm1, 1);
                     rowmax = max(rowmax, abs(real(Wkwm1[jmax]))+abs(imag(Wkwm1[jmax])));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                     kp = k;
                  else if (abs(real(Wkwm1[imax]))+abs(imag(Wkwm1[imax])) >= alpha*rowmax)
                  {
                     kp = imax;
                     LATL::COPY(k+1, Wkwm1, 1, Wkw, 1);
                  }
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               
               kk = k - kstep +1;
               kkw = nb + kk -n;
               Akk = A + ldA*kk;
               Wkkw = Work + ldWork*kkw;
               if (kp != kk)
               {
                  Ak[kp] = Ak[kk];
                  LATL::COPY(k-kp-1, Akk+kp+1, 1, A + ldA*(kp+1)+ kp, ldA);
                  LATL::COPY(kp+1, Akk, 1, A+ldA*kp, 1);
                  LATL::SWAP(n-kk, Akk+kk, ldA, Akk+kp, ldA);
                  LATL::SWAP(n-kk, Wkkw+kk, ldWork, Wkkw+kp, ldWork);
               }
               
               if (kstep == 1)
               {
                  LATL::COPY(k+1, Wkw, 1, Ak, 1);
                  r1 = onec/Ak[k];
                  LATL::SCAL(k, r1, Ak, 1);
               }
               else
               {
                  if (k > 1)
                  {
                     d21 = Wkw[k-1];
                     d11 = Wkw[k]/d21;
                     d22 = Wkwm1[k-1]/d21;
                     d21 = (onec/(d11*d22-onec))/d21;
                     for (int_t j = 0; j <= k-2; ++j)
                     {
                        Akm1[j] = d21*(d11*Wkwm1[j]-Wkw[j]);
                        Ak[j] = d21*(d22*Wkw[j] - Wkwm1[j]);
                     }
                  }
                  
                  Akm1[k-1] = Wkwm1[k-1];
                  Ak[k-1] = Wkw[k-1];
                  Ak[k] = Wkw[k];
               }
               
            }
            if (kstep == 1)
            {
               ipiv[k] = kp;
               bsdv[k] = 0;
            }
            else
            {
               ipiv[k] = kp;
               ipiv[k-1] = kp;
               bsdv[k] = 1;
               bsdv[k-1] = 1;
            }
            
            k -= kstep;
         }
         
         for (int_t j = ((k)/nb)*nb; j >= 0; j -= nb)
         {
            jb = min(nb, k-j+1);
            for (int_t jj = j; jj < j+jb; ++jj)
            {
               LATL::GEMV('N', jj-j+1, n-k-1,-onec, Ak+ldA+j, ldA, Wkw+ldWork+jj, ldWork, onec, A+jj*ldA+j, 1);
            }
            LATL::GEMM('N', 'T', j, jb, n-k-1, -onec, Ak+ldA, ldA, Wkw+ldWork+j, ldWork, onec, A+j*ldA, ldA);
         }
         
         kp1 = k+1;
         while (kp1 < n)
         {
            jtemp = kp1;
            jp = ipiv[kp1];
            if (bsdv[kp1] == 1)
            {
               ++kp1;
            }
            ++kp1;
            if (jp != jtemp && kp1 < n)
            {
               LATL::SWAP(n-kp1, A+ldA*kp1+jp, ldA, A+ldA*kp1+jtemp, ldA);
            }
         }
         kb = n-k-1;
         
      }
      else
      {
         k = 0;
         complex<real_t> * Wk, *Wkp1, *Akp1;
         int_t km1;
         while (k < n)
         {
            if (k >= nb-1 && nb < n)
               break;
            
            Ak = A+ldA*k;
            Akp1 = Ak+ldA;
            Wk = Work+ldWork*k;
            Wkp1 = Wk+ldWork;
            LATL::COPY(n-k, Ak+k, 1, Wk+k, 1);
            LATL::GEMV('N', n-k, k, -onec, A+k, ldA, Work+k, ldWork, onec, Wk+k, 1);
            
            kstep = 1;
            
            absakk = abs(real(Wk[k]))+abs(imag(Wk[k]));
            
            if (k < n-1)
            {
               imax = k+1 + LATL::IMAX(n-k-1, Wk+k+1, 1);
               colmax = abs(real(Wk[imax]))+abs(imag(Wk[imax]));
            }
            else
            {
               colmax = zero;
            }
            
            if (max(absakk, colmax) == zero)
            {
               if (info == 0)
                  info = k+1;
               kp = k;
            }
            else
            {
               if (absakk >= alpha*colmax)
                  kp = k;
               else
               {
                  Aimax = A+ldA*imax;
                  LATL::COPY(imax-k, Ak+imax, ldA, Wkp1+k, 1);
                  LATL::COPY(n-imax, Aimax+imax, 1, Wkp1+imax, 1);
                  LATL::GEMV('N', n-k, k, -onec, A+k, ldA, Work+imax, ldWork, onec, Wkp1+k, 1);
                  
                  jmax = k + LATL::IMAX(imax-k, Wkp1+k, 1);
                  rowmax = abs(real(Wkp1[jmax]))+abs(imag(Wkp1[jmax]));
                  if (imax < n-1)
                  {
                     jmax = imax+1+LATL::IMAX(n-imax-1, Wkp1+imax+1, 1);
                     rowmax = max(rowmax, abs(real(Wkp1[jmax]))+abs(imag(Wkp1[jmax])));
                  }
                  
                  if (absakk >= alpha*colmax*(colmax/rowmax))
                  {
                     kp = k;
                  }
                  else if ( abs(real(Wkp1[imax]))+abs(imag(Wkp1[imax])) >= alpha*rowmax)
                  {
                     kp = imax;
                     LATL::COPY(n-k, Wkp1+k, 1, Wk+k, 1);
                  }
                  else
                  {
                     kp = imax;
                     kstep = 2;
                  }
               }
               kk = k+kstep-1;
               Akk = A+ldA*kk;
               
               if (kp != kk)
               {
                  Ak[kp] = Ak[kk];
                  LATL::COPY(kp-k-1, Akk+k+1, 1, Akp1+kp, ldA);
                  LATL::COPY(n-kp, Akk+kp, 1, A+ldA*kp+kp, 1);
                  LATL::SWAP(kk+1, A+kk, ldA, A+kp, ldA);
                  LATL::SWAP(kk+1, Work+kk, ldWork, Work+kp, ldWork);
               }
               
               if (kstep == 1)
               {
                  LATL::COPY(n-k, Wk+k, 1, Ak+k, 1);
                  if (k < n-1)
                  {
                     r1 = onec/Ak[k];
                     LATL::SCAL(n-k-1, r1, Ak+k+1, 1);
                  }
               }
               else
               {
                  if ( k < n-2)
                  {
                     d21 = Wk[k+1];
                     d11 = Wkp1[k+1]/d21;
                     d22 = Wk[k]/d21;
                     d21 = (onec/(d11*d22-onec))/d21;
                     for (int_t j = k+2; j < n; ++j)
                     {
                        Ak[j] = d21*(d11*Wk[j] - Wkp1[j]);
                        Akp1[j] = d21*(d22*Wkp1[j]-Wk[j]);
                     }
                  }
                  
                  Ak[k] = Wk[k];
                  Ak[k+1] = Wk[k+1];
                  Akp1[k+1] = Wkp1[k+1];
               }
            }
            
            if (kstep == 1)
            {
               ipiv[k] = kp;
               bsdv[k] = 0;
            }
            else
            {
               ipiv[k] = kp;
               ipiv[k+1] = kp;
               bsdv[k] = 1;
               bsdv[k+1] = 1;
            }
            
            k += kstep;
         }
         for (int_t j = k; j < n; j += nb)
         {
            jb = min(nb, n-j);
            for (int_t jj = j; jj < j+jb; ++jj)
            {
               LATL::GEMV('N', j+jb-jj, k, -onec, A+jj, ldA, Work+jj, ldWork, onec, A+jj*ldA+jj, 1);
            }
            
            if (j+jb < n)
            {
               LATL::GEMM('N', 'T', n-j-jb, jb, k, -onec, A+j+jb, ldA, Work+j, ldWork, onec, A+j*ldA+j+jb, ldA);
            }
         }
         
         km1 = k-1;
         while (km1 >= 0)
         {
            jtemp = km1;
            jp = ipiv[km1];
            if (bsdv[km1] == 1)
            {
               --km1;
            }
            --km1;
            if (jp != jtemp && km1 >= 0)
            {
               LATL::SWAP(km1+1, A+jp, ldA, A+jtemp, ldA);
            }
         }
         kb = k;
      }
      
      return info;
   }
}


#endif
