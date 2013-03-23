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

#include <limits>
#include "imax.h"
#include "gemv.h"
#include "scal.h"
#include "swap.h"
#include "lacgv.h"
#include "syrk.h"
#include "herk.h"
#include "latl.h"

namespace LATL
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
   /// @ingroup TRF

   template< typename real_t>
   int PSTRF(const char uplo, const int_t n, real_t * const A, const int_t ldA, int_t * const PIV, int_t &rank, const real_t tol)
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

      for (int_t i = 0; i < n; ++i)
         PIV[i] = i;

      const real_t one(1.0);
      const real_t zero(0.0);
      real_t ajj, stop, temp, *Ai = A+ldA, *Aj = A, *Apvt, *Apvtp1;
      int_t itemp, pvt = 0;

      ajj = A[0];
      for (int_t i = 1; i < n; ++i)
      {
         if (Ai[i] > ajj)
         {
            pvt = i;
            ajj = Ai[i];
         }
         Ai += ldA;
      }
      if (ajj == zero || std::isnan(ajj))
      {
         rank = 0;
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

      real_t * Work = new real_t[2*n];
      for (int_t i = 0; i < n; ++i)
      {
         Work[i] = 0;
      }

      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = j; i < n; ++i)
            {
               Ai = A+ldA*i;
               if (j > 0)
                  Work[i] += Ai[j-1]*Ai[j-1];
               Work[n+i] = Ai[i] - Work[i];
            }

            if ( j > 0)
            {
               pvt = j;
               for (int_t i =j+1; i < n; ++i)
               {
                  if (Work[n+i] > Work[n+pvt])
                     pvt = i;
               }
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
               Apvtp1 = Apvt+ldA;
               Apvt[pvt] = Aj[j];
               LATL::SWAP(j, Aj, 1, Apvt, 1);
               if (pvt < n-1)
               {
                  LATL::SWAP(n-pvt-1, Apvtp1+j, ldA, Apvtp1+pvt, ldA);
               }
               LATL::SWAP(pvt-j-1, Aj+ldA+j, ldA, Apvt+j+1, 1);

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
               LATL::GEMV('T', j, n-j-1, -one, Aj+ldA, ldA, Aj, 1, one, Aj+ldA+j, ldA);
               LATL::SCAL(n-j-1, one/ajj, Aj+ldA+j, ldA);
            }
            Aj += ldA;
         }
      }
      else
      {
         real_t *Ajm1;
         for (int_t j = 0; j < n; ++j)
         {
            Ajm1 = Aj - ldA;
            for (int_t i = j; i < n; ++i)
            {
               Ai = A + ldA*i;
               if (j > 0)
               {
                  Work[i] += Ajm1[i]*Ajm1[i];
               }
               Work[n+i] = Ai[i] - Work[i];
            }

            if ( j > 0)
            {
               pvt = j;
               for (int_t i =j+1; i < n; ++i)
               {
                  if (Work[n+i] > Work[n+pvt])
                     pvt = i;
               }
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
               LATL::SWAP(j, A+j, ldA, A+pvt, ldA);
               if (pvt < n-1)
               {
                  LATL::SWAP(n-pvt-1, Aj+pvt+1, 1, Apvt+pvt+1, 1);
               }
               LATL::SWAP(pvt-j-1, Aj+j+1, 1, Aj+ldA+pvt, ldA);

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
               LATL::GEMV('N', n-j-1, j, -one, A+j+1, ldA, A+j, ldA, one, Aj+j+1, 1);
               LATL::SCAL(n-j-1, one/ajj, Aj+j+1, 1);
            }
            Aj += ldA;
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
   /// @param A Complex array size ldA-by-n.  On entry, the Hermitian matrix A.  On exit, if the return value is 0, the factor U or L from the Cholesky factorization.
   /// @param ldA Column length of the matrix A.  ldA >= n
   /// @param PIV Integer array size n.  The nonzero entries of PIV indicate P(PIV[k], k) = 1
   /// @param rank Integer.  On exit, indicates the number of steps the algorithm completed.
   /// @param tol User defined tolerance.  If tol < 0, the n*U*max(Ak[k]) will be used.  The algorithm terminates at the k-1st step if the pivot <= tol.
   /// @ingroup TRF

   template< typename real_t> int_t PSTRF(const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA, int_t * const PIV, int_t &rank, const real_t tol)
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

      for (int_t i = 0; i < n; ++i)
         PIV[i] = i;

      const real_t one(1.0);
      const complex<real_t> onec(1.0);
      const real_t zero(0.0);
      real_t ajj, stop, temp;
      complex<real_t> *Ai = A, *Aj = A, *Apvt, *Apvtp1, ztemp;
      int_t itemp, pvt = 0;
      real_t * Work = new real_t[2*n];

      for (int_t i = 0; i < n; ++i)
      {
         Work[i] = real(Ai[i]);
         Ai += ldA;
      }
      pvt = 0;
      for (int_t i =1; i < n; ++i)
      {
         if (Work[i] > Work[pvt])
            pvt = i;
      }
      Apvt = A + ldA*pvt;
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

      for (int_t i = 0; i < n; ++i)
      {
         Work[i] = 0;
      }

      if (uplo == 'U' || uplo == 'u')
      {
         for (int_t j = 0; j < n; ++j)
         {
            for (int_t i = j; i < n; ++i)
            {
               Ai = A+ldA*i;
               if (j > 0)
                  Work[i] += real((conj(Ai[j-1])*Ai[j-1]));//
               Work[n+i] = real(Ai[i]) - Work[i];//
            }

            if ( j > 0)
            {
               pvt = j;
               for (int_t i =j+1; i < n; ++i)
               {
                  if (Work[n+i] > Work[n+pvt])
                     pvt = i;
               }
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
               Apvtp1 = Apvt+ldA;
               Apvt[pvt] = Aj[j];
               LATL::SWAP(j, Aj, 1, Apvt, 1);
               if (pvt < n-1)
               {
                  LATL::SWAP(n-pvt-1, Apvtp1+j, ldA, Apvtp1+pvt, ldA);
               }
               for (int_t i = j+1; i < pvt; ++i)//
               {
                  Ai = A+ldA*i;
                  ztemp = conj(Ai[j]);
                  Ai[j] = conj(Apvt[i]);
                  Apvt[i] = ztemp;
               }
               Apvt[j] = conj(Apvt[j]);

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
               LATL::LACGV(j, Aj, 1);//
               LATL::GEMV('T', j, n-j-1, -onec, Aj+ldA, ldA, Aj, 1, onec, Aj+ldA+j, ldA);
               LATL::LACGV(j, Aj, 1);//
               LATL::SCAL(n-j-1, one/ajj, Aj+ldA+j, ldA);
            }
            Aj += ldA;
         }
      }
      else
      {
         complex<real_t> *Ajm1;
         for (int_t j = 0; j < n; ++j)
         {
            Ajm1 = Aj - ldA;
            for (int_t i = j; i < n; ++i)
            {
               Ai = A + ldA*i;
               if (j > 0)
               {
                  Work[i] += real(conj(Ajm1[i])*Ajm1[i]);
               }
               Work[n+i] = real(Ai[i]) - Work[i];

            }

            if ( j > 0)
            {
               pvt = j;
               for (int_t i =j+1; i < n; ++i)
               {
                  if (Work[n+i] > Work[n+pvt])
                     pvt = i;
               }
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
               LATL::SWAP(j, A+j, ldA, A+pvt, ldA);
               if (pvt < n-1)
               {
                  LATL::SWAP(n-pvt-1, Aj+pvt+1, 1, Apvt+pvt+1, 1);
               }
               for (int_t i = j+1; i < pvt; ++i)
               {
                  Ai = A + ldA*i;
                  ztemp = conj(Aj[i]);
                  Aj[i] = conj(Ai[pvt]);
                  Ai[pvt] = ztemp;
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
               LATL::LACGV(j, A+j, ldA);
               LATL::GEMV('N', n-j-1, j, -onec, A+j+1, ldA, A+j, ldA, onec, Aj+j+1, 1);
               LATL::LACGV(j, A+j, ldA);
               LATL::SCAL(n-j-1, one/ajj, Aj+j+1, 1);
            }
            Aj += ldA;
         }
      }

      rank = n;
      delete [] Work;
      return 0;
   }

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
   /// @param nb Block size.
   
   template< typename real_t>
   int PSTRF(const char uplo, const int_t n, real_t * const A, const int_t ldA, int_t * PIV, int_t &rank, const real_t tol, const int_t nb)
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
         return LATL::PSTRF(uplo, n, A, ldA, PIV, rank, tol);
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
                  itemp = LATL::IMAX(n-j , Work+n+j, 1);
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
                  LATL::SWAP(j, Aj, 1, Apvt, 1);
                  if (pvt < n-1)
                  {
                     LATL::SWAP(n-pvt-1, Apvt+ldA+j, ldA, Apvt+ldA+pvt, ldA);
                  }
                  LATL::SWAP(pvt-j-1, Ajp1+j, ldA, Apvt+j+1, 1);
                  
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
                  LATL::GEMV('T', j-k, n-j-1, -one, Ajp1+k, ldA, Aj+k, 1, one, Ajp1+j, ldA);
                  LATL::SCAL(n-j-1, one/ajj, Ajp1+j, ldA);
               }
               Aj += ldA;
               Ajp1 += ldA;
               
            }
            
            itemp = k+jb;
            Aj = A + ldA*itemp;
            if (itemp < n)
            {
               LATL::SYRK('U', 'T', n-itemp, jb, -one, Aj+k, ldA, one, Aj+itemp, ldA);
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
                  itemp = LATL::IMAX(n-j , Work+n+j, 1);
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
                  LATL::SWAP(j, A+j, ldA, A+pvt, ldA);
                  if (pvt < n-1)
                  {
                     LATL::SWAP(n-pvt-1, Aj+pvt+1, 1, Apvt+pvt+1, 1);
                  }
                  LATL::SWAP(pvt-j-1, Aj+j+1, 1, Aj+ldA+pvt, ldA);
                  
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
                  LATL::GEMV('N', n-j-1, j-k, -one, Ak+j+1, ldA, Ak+j, ldA, one, Aj+j+1, 1);
                  LATL::SCAL(n-j-1, one/ajj, Aj+j+1, 1);
               }
               
               Aj += ldA;
               Ajm1 += ldA;
            }
            
            itemp = k+jb;
            Aj = A + ldA*itemp;
            if (itemp < n)
            {
               LATL::SYRK('L', 'N', n-itemp, jb, -one, Ak+itemp, ldA, one, Aj+itemp, ldA);
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
   /// @param nb Block size.

   template< typename real_t>
   int PSTRF(const char uplo, const int_t n, complex<real_t> * const A, const int_t ldA, int_t * PIV, int_t &rank, const real_t tol, const int_t nb)
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
         return LATL::PSTRF(uplo, n, A, ldA, PIV, rank, tol);
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
      
      pvt = LATL::IMAX(n, Work, 1);
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
                  itemp = LATL::IMAX(n-j , Work+n+j, 1);
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
                  LATL::SWAP(j, Aj, 1, Apvt, 1);
                  if (pvt < n-1)
                  {
                     LATL::SWAP(n-pvt-1, Apvt+ldA+j, ldA, Apvt+ldA+pvt, ldA);
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
                  LATL::LACGV(j, Aj, 1);
                  LATL::GEMV('T', j-k, n-j-1, -onec, Ajp1+k, ldA, Aj+k, 1, onec, Ajp1+j, ldA);
                  LATL::LACGV(j, Aj, 1);
                  LATL::SCAL(n-j-1, one/ajj, Ajp1+j, ldA);
               }
               Aj += ldA;
               Ajp1 += ldA;
            }
            
            itemp = k+jb;
            Aj = A + ldA*itemp;
            if (itemp < n)
            {
               LATL::HERK('U', 'C', n-itemp, jb, -one, Aj+k, ldA, one, Aj+itemp, ldA);
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
                  itemp = LATL::IMAX(n-j , Work+n+j, 1);
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
                  LATL::SWAP(j, A+j, ldA, A+pvt, ldA);
                  if (pvt < n-1)
                  {
                     LATL::SWAP(n-pvt-1, Aj+pvt+1, 1, Apvt+pvt+1, 1);
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
                  LATL::LACGV(j, A+j, ldA);
                  LATL::GEMV('N', n-j-1, j-k, -onec, Ak+j+1, ldA, Ak+j, ldA, onec, Aj+j+1, 1);
                  LATL::LACGV(j, A+j, ldA);
                  LATL::SCAL(n-j-1, one/ajj, Aj+j+1, 1);
               }
               
               Aj += ldA;
               Ajm1 += ldA;
            }
            
            itemp = k+jb;
            Aj = A + ldA*itemp;
            if (itemp < n)
            {
               LATL::HERK('L', 'N', n-itemp, jb, -one, Ak+itemp, ldA, one, Aj+itemp, ldA);
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
