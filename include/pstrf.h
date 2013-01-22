//
//  pstrf.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 8/1/12.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _pstrf_h
#define _pstrf_h

/// @file pstrf.h Computes the Cholesky factorization with complete pivoting of a positive semidefinite matrix.

#include "pstf2.h"
#include "gemv.h"
#include "scal.h"
#include "swap.h"
#include "syrk.h"
#include "imax.h"
#include "lacgv.h"
#include "herk.h"
#include <limits>
#include "latl.h"

namespace latl
{
   /// @brief Computes the Cholesky factorization with complete pivoting of a real symmetric positive semidefinite matrix A.
   ///
   /// The factorization has the form
   ///
   ///     P' * A * P = U' * U     if uplo = 'U'
   ///     P' * A * P = L * L'     if uplo = 'L'
   ///
   /// where U is an upper triangular matrix and L is lower triangular, and P is stored as vector PIV.
   ///
   /// This algorithm does not attempt to check that A is positive semidefinite.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return 1 if the matrix A is rank deficient, in which case the computed rank is returned in the argument rank, or is indefinite.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetrix matrix A is stored.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Real array size ldA-by-n.  On entry, the symmetric matrix A.  On exit, if the return value is 0, the factor U or L from the Cholesky factorization.
   /// @param ldA Column length of the matrix A.  ldA >= n
   /// @param PIV Integer array size n.  The nonzero entries of PIV indicate P(PIV[k], k) = 1
   /// @param rank Integer.  On exit, indicates the number of steps the algorithm completed.
   /// @param tol User defined tolerance.  If tol < 0, the n*U*max(Ak[k]) will be used.  The algorithm terminates at the k-1st step if the pivot <= tol.
   /// @param nb Block size, optional.  Default value is 64.
   /// @ingroup TRF
   
   template< typename real_t>
   int pstrf(const char uplo, const int_t n, real_t * const A, const int_t ldA, int_t * PIV, int_t &rank, const real_t tol, const int_t nb = 64)
   {
      using std::numeric_limits;

      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      if (nb <= 1 || nb >= n)
      {
         return latl::pstf2(uplo, n, A, ldA, PIV, rank, tol);
      }
      for (int_t i = 0; i < n; ++i)
      {
         PIV[i] = i;
      }
      
      const real_t one(1.0);
      const real_t zero(0.0);
      real_t ajj, stop, temp, *Apvt = A, *Ai = A + ldA, *Aj, *Ajp1, *Ajm1, *Ak;
      int_t pvt = 0, itemp, jb;
      
      ajj = Apvt[pvt];
      for (int_t i = 1; i < n; ++i)
      {
         if (Ai[i] > ajj)
         {
            pvt = i;
            Apvt = A + ldA*pvt;
            ajj = Apvt[pvt];
         }
         Ai += ldA;
      }
      
      if (ajj == zero || std::isnan(ajj))
      {
         rank = 0;
         return 1;
      }
      
      real_t *Work = new real_t[2*n];
      if (tol < zero)
      {
         stop = n * numeric_limits<real_t>::epsilon()*ajj;
      }
      else
      {
         stop = tol;
      }
      
      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t k = 0; k < n; k += nb)
         {
            jb = std::min(nb, n-k);
            
            for (int_t i = k; i < n; ++i)
            {
               Work[i] = 0;
            }
            
            Aj = A + ldA*k;
            Ajp1 = Aj+ldA;
            for (int_t j = k; j < k+jb; ++j)
            {
               
               for (int_t i = j; i < n; ++i)
               {
                  Ai = A + ldA*i;
                  if ( j > k)
                  {
                     Work[i] += Ai[j-1]*Ai[j-1];
                  }
                  Work[n+i] = Ai[i] - Work[i];
               }
               
               if (j > 0)
               {
                  itemp = latl::imax(n-j , Work+n+j, 1);
                  pvt = itemp + j;
                  ajj = Work[n+pvt];
                  if (ajj <= stop || std::isnan(ajj))
                  {
                     Aj[j] = ajj;
                     rank = j;
                     delete [] Work;
                     return 1;
                  }
               }
               
               if (j != pvt)
               {
                  Apvt = A + ldA*pvt;
                  Apvt[pvt] = Aj[j];
                  latl::swap(j, Aj, 1, Apvt, 1);
                  if (pvt < n-1)
                  {
                     latl::swap(n-pvt-1, Apvt+ldA+j, ldA, Apvt+ldA+pvt, ldA);
                  }
                  latl::swap(pvt-j-1, Ajp1+j, ldA, Apvt+j+1, 1);
                  
                  temp = Work[j];
                  Work[j] = Work[pvt];
                  Work[pvt] = temp;
                  itemp = PIV[pvt];
                  PIV[pvt] = PIV[j];
                  PIV[j] = itemp;
               }
               
               Aj[j] = ajj = std::sqrt(ajj);
               
               if (j < n-1)
               {
                  latl::gemv('T', j-k, n-j-1, -one, Ajp1+k, ldA, Aj+k, 1, one, Ajp1+j, ldA);
                  latl::scal(n-j-1, one/ajj, Ajp1+j, ldA);
               }
               Aj += ldA;
               Ajp1 += ldA;
               
            }
            
            itemp = k+jb;
            Aj = A + ldA*itemp;
            if (itemp < n)
            {
               latl::syrk('U', 'T', n-itemp, jb, -one, Aj+k, ldA, one, Aj+itemp, ldA);
            }
         }
      }
      else
      {
         Ak = A;
         for (int_t k = 0; k < n; k += nb)
         {
            jb = std::min(nb, n-k);
            
            for (int_t i = k; i < n; ++i)
            {
               Work[i] = 0;
            }
            
            Aj = A +ldA*k;
            Ajm1 = Aj - ldA;
            for (int_t j = k; j < k+jb; ++j)
            {
               Ai = A+ldA*j;
               for (int_t i = j; i < n; ++i)
               {
                  if (j > k)
                  {
                     Work[i] += Ajm1[i]*Ajm1[i];
                  }
                  Work[n+i] = Ai[i] - Work[i];
                  Ai += ldA;
               }
               
               if (j > 0)
               {
                  itemp = latl::imax(n-j , Work+n+j, 1);
                  pvt = itemp + j;
                  ajj = Work[n+pvt];
                  if (ajj <= stop || std::isnan(ajj))
                  {
                     Aj[j] = ajj;
                     rank = j;
                     delete [] Work;
                     return 1;
                  }
               }
               
               if (j != pvt)
               {
                  Apvt = A+ldA*pvt;
                  Apvt[pvt] = Aj[j];
                  latl::swap(j, A+j, ldA, A+pvt, ldA);
                  if (pvt < n-1)
                  {
                     latl::swap(n-pvt-1, Aj+pvt+1, 1, Apvt+pvt+1, 1);
                  }
                  latl::swap(pvt-j-1, Aj+j+1, 1, Aj+ldA+pvt, ldA);
                  
                  temp = Work[j];
                  Work[j] = Work[pvt];
                  Work[pvt] = temp;
                  itemp = PIV[pvt];
                  PIV[pvt] = PIV[j];
                  PIV[j] = itemp;
               }
               
               Aj[j] = ajj = std::sqrt(ajj);
               
               if (j < n-1)
               {
                  latl::gemv('N', n-j-1, j-k, -one, Ak+j+1, ldA, Ak+j, ldA, one, Aj+j+1, 1);
                  latl::scal(n-j-1, one/ajj, Aj+j+1, 1);
               }
               
               Aj += ldA;
               Ajm1 += ldA;
            }
            
            itemp = k+jb;
            Aj = A + ldA*itemp;
            if (itemp < n)
            {
               latl::syrk('L', 'N', n-itemp, jb, -one, Ak+itemp, ldA, one, Aj+itemp, ldA);
            }
            
            Ak += ldA*nb;
         }
      }
      
      rank = n;
      delete [] Work;
      return 0;
   }
   
   /// @brief Computes the Cholesky factorization with complete pivoting of a complex Hermitian positive semidefinite matrix A.
   ///
   /// The factorization has the form
   ///
   ///     P' * A * P = U' * U     if uplo = 'U'
   ///     P' * A * P = L * L'     if uplo = 'L'
   ///
   /// where U is an upper triangular matrix and L is lower triangular, and P is stored as vector PIV.
   ///
   /// This algorithm does not attempt to check that A is positive semidefinite.
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @return 1 if the matrix A is rank deficient, in which case the computed rank is returned in the argument rank, or is indefinite.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether the upper or lower triangular part of the symmetrix matrix A is stored.
   /// @param n Order of the matrix A.  n >= 0
   /// @param A Real array size ldA-by-n.  On entry, the symmetric matrix A.  On exit, if the return value is 0, the factor U or L from the Cholesky factorization.
   /// @param ldA Column length of the matrix A.  ldA >= n
   /// @param PIV Integer array size n.  The nonzero entries of PIV indicate P(PIV[k], k) = 1
   /// @param rank Integer.  On exit, indicates the number of steps the algorithm completed.
   /// @param tol User defined tolerance.  If tol < 0, the n*U*max(Ak[k]) will be used.  The algorithm terminates at the k-1st step if the pivot <= tol.
   /// @param nb Block size, optional.  Default value is 64.
   /// @ingroup TRF
   
   template< typename real_t>
   int pstrf(const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA, int_t * PIV, int_t &rank, const real_t tol, const int_t nb = 64)
   {
      using std::numeric_limits;

      if (uplo != 'U' && uplo != 'L' && uplo != 'u' && uplo != 'l')
         return -1;
      if (n < 0)
         return -2;
      if (ldA < n)
         return -4;
      
      if (n == 0)
         return 0;
      
      if (nb <= 1 || nb >= n)
      {
         return latl::pstf2(uplo, n, A, ldA, PIV, rank, tol);
      }
      for (int_t i = 0; i < n; ++i)
      {
         PIV[i] = i;
      }
      
      const real_t one(1.0);
      const complex<real_t> onec(1.0);
      const real_t zero(0.0);
      real_t ajj, stop, temp;
      complex<real_t> *Apvt = A, *Ai = A, *Aj, *Ajp1, *Ajm1, *Ak, ztemp;
      int_t pvt, itemp, jb;
      real_t *Work = new real_t[2*n];
      
      for (int_t i = 0; i < n; ++i)
      {
         Work[i] = real(Ai[i]);
         Ai += ldA;
      }
      
      pvt = latl::imax(n, Work, 1);
      Apvt += ldA*pvt;
      ajj = real(Apvt[pvt]);
      
      if (ajj == zero || std::isnan(ajj))
      {
         rank = 0;
         delete [] Work;
         return 1;
      }
      
      if (tol < zero)
      {
         stop = n * numeric_limits<real_t>::epsilon()*ajj;
      }
      else
      {
         stop = tol;
      }
      
      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t k = 0; k < n; k += nb)
         {
            jb = std::min(nb, n-k);
            
            for (int_t i = k; i < n; ++i)
            {
               Work[i] = 0;
            }
            
            Aj = A + ldA*k;
            Ajp1 = Aj+ldA;
            for (int_t j = k; j < k+jb; ++j)
            {
               
               for (int_t i = j; i < n; ++i)
               {
                  Ai = A + ldA*i;
                  if ( j > k)
                  {
                     Work[i] += real(conj(Ai[j-1])*Ai[j-1]);
                  }
                  Work[n+i] = real(Ai[i]) - Work[i];
               }
               
               if (j > 0)
               {
                  itemp = latl::imax(n-j , Work+n+j, 1);
                  pvt = itemp + j;
                  ajj = Work[n+pvt];
                  if (ajj <= stop || std::isnan(ajj))
                  {
                     Aj[j] = ajj;
                     rank = j;
                     delete [] Work;
                     return 1;
                  }
               }
               
               if (j != pvt)
               {
                  Apvt = A + ldA*pvt;
                  Apvt[pvt] = Aj[j];
                  latl::swap(j, Aj, 1, Apvt, 1);
                  if (pvt < n-1)
                  {
                     latl::swap(n-pvt-1, Apvt+ldA+j, ldA, Apvt+ldA+pvt, ldA);
                  }
                  Ai = A+ldA*(j+1);
                  for (int_t i = j+1; i < pvt; ++i)
                  {
                     ztemp = conj(Ai[j]);
                     Ai[j] = conj(Apvt[i]);
                     Apvt[i] = ztemp;
                     Ai += ldA;
                  }
                  
                  temp = Work[j];
                  Work[j] = Work[pvt];
                  Work[pvt] = temp;
                  itemp = PIV[pvt];
                  PIV[pvt] = PIV[j];
                  PIV[j] = itemp;
               }
               
               Aj[j] = ajj = std::sqrt(ajj);
               
               if (j < n-1)
               {
                  latl::lacgv(j, Aj, 1);
                  latl::gemv('T', j-k, n-j-1, -onec, Ajp1+k, ldA, Aj+k, 1, onec, Ajp1+j, ldA);
                  latl::lacgv(j, Aj, 1);
                  latl::scal(n-j-1, one/ajj, Ajp1+j, ldA);
               }
               Aj += ldA;
               Ajp1 += ldA;
            }
            
            itemp = k+jb;
            Aj = A + ldA*itemp;
            if (itemp < n)
            {
               latl::herk('U', 'C', n-itemp, jb, -one, Aj+k, ldA, one, Aj+itemp, ldA);
            }
         }
      }
      else
      {
         Ak = A;
         for (int_t k = 0; k < n; k += nb)
         {
            jb = std::min(nb, n-k);
            
            for (int_t i = k; i < n; ++i)
            {
               Work[i] = 0;
            }
            
            Aj = A +ldA*k;
            Ajm1 = Aj - ldA;
            for (int_t j = k; j < k+jb; ++j)
            {
               Ai = A+ldA*j;
               for (int_t i = j; i < n; ++i)
               {
                  if (j > k)
                  {
                     Work[i] += real(conj(Ajm1[i])*Ajm1[i]);
                  }
                  Work[n+i] = real(Ai[i]) - Work[i];
                  Ai += ldA;
               }
               
               if (j > 0)
               {
                  itemp = latl::imax(n-j , Work+n+j, 1);
                  pvt = itemp + j;
                  ajj = Work[n+pvt];
                  if (ajj <= stop || std::isnan(ajj))
                  {
                     Aj[j] = ajj;
                     rank = j;
                     delete [] Work;
                     return 1;
                  }
               }
               
               if (j != pvt)
               {
                  Apvt = A+ldA*pvt;
                  Apvt[pvt] = Aj[j];
                  latl::swap(j, A+j, ldA, A+pvt, ldA);
                  if (pvt < n-1)
                  {
                     latl::swap(n-pvt-1, Aj+pvt+1, 1, Apvt+pvt+1, 1);
                  }
                  Ai = A + ldA*(j+1);
                  for (int_t i = j+1; i < pvt; ++i)
                  {
                     ztemp = Aj[i];
                     Aj[i] = conj(Ai[pvt]);
                     Ai[pvt] = ztemp;
                     Ai += ldA;
                  }
                  Aj[pvt] = conj(Aj[pvt]);
                  
                  temp = Work[j];
                  Work[j] = Work[pvt];
                  Work[pvt] = temp;
                  itemp = PIV[pvt];
                  PIV[pvt] = PIV[j];
                  PIV[j] = itemp;
               }
               
               Aj[j] = ajj = std::sqrt(ajj);
               
               if (j < n-1)
               {
                  latl::lacgv(j, A+j, ldA);
                  latl::gemv('N', n-j-1, j-k, -onec, Ak+j+1, ldA, Ak+j, ldA, onec, Aj+j+1, 1);
                  latl::lacgv(j, A+j, ldA);
                  latl::scal(n-j-1, one/ajj, Aj+j+1, 1);
               }
               
               Aj += ldA;
               Ajm1 += ldA;
            }
            
            itemp = k+jb;
            Aj = A + ldA*itemp;
            if (itemp < n)
            {
               latl::herk('L', 'N', n-itemp, jb, -one, Ak+itemp, ldA, one, Aj+itemp, ldA);
            }
            
            Ak += ldA*nb;
         }
      }
      
      rank = n;
      delete [] Work;
      return 0;
   }
   
}


#endif
