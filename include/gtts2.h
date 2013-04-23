//
//  gtts2.h
//  Linear Algebra Template Library
//
//  Created by Stephanie Patterson on 3/22/13.
//  Copyright (c) 2012 University of Colorado Denver. All rights reserved.
//

#ifndef _gtts2_h
#define _gtts2_h

/// @file gtts2.h Solves a system of linear equations A * X = B.

#include "latl.h"

namespace LATL
{
   /// @brief Solves a real system of equations of the form
   ///
   ///      A * X = B      or    A' * X = B
   ///
   /// using the LU factorization of the tridiagonal matrix A computed by GTTRF.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param itrans Specifies the form of the system of equations.
   ///
   ///      0 = no transpose
   ///      1 = transpose
   ///
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param DL Real array, size n-1.  On entry, should contain the (n-1) multipliers that define the matrix L from the LU factorization of A.
   /// @param D Real array, size n.  On entry, should contain the n diagonal elements of the upper triangular matrix U from the LU factorization of A.
   /// @param DU Real array, size n-1.  On entry, the (n-1) elements of the first superdiagonal of U.
   /// @param DU2 Real array, size n-2.  On entry, the (n-2) elements of the second superdiagonal of U.
   /// @param ipiv Integer array, size n.  The pivot indices from the factorization LU.
   /// @param B Real matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B.
   /// @ingroup COMP
   
   template< typename real_t>
   int_t GTTS2(const int_t itrans, const int_t n, const int_t nrhs, real_t * const DL, real_t * const D, real_t * const DU, real_t * const DU2, int_t * const ipiv, real_t * const B, const int_t ldB)
   {
      if (itrans != 0 && itrans != 1 && itrans != 2)
         return -1;
      if (n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldB < n)
         return -10;
      if (n == 0 || nrhs == 0)
         return 0;
      if (itrans == 0)
      {
         if (nrhs == 1)
         {
            for (int_t i = 0; i < n-1; ++i)
            {
               int_t ip = ipiv[i];
               real_t temp = B[2*i-ip+1]-DL[i]*B[ip];
               B[i] = B[ip];
               B[i+1] = temp;
            }
            B[n-1] = B[n-1]/D[n-1];
            if (n > 1)
            {
               B[n-2] = (B[n-2]-DU[n-2]*B[n-1])/D[n-2];
            }
            for (int_t i = n-3; i >= 0; --i)
            {
               B[i] = (B[i]-DU[i]*B[i+1]-DU2[i]*B[i+2])/D[i];
            }
         }
         else
         {
            real_t * Bj = B;
            for (int_t j = 0; j < nrhs; ++j)
            {
               for (int_t i = 0; i < n-1; ++i)
               {
                  if (ipiv[i] == i)
                     Bj[i+1] -= DL[i]*Bj[i];
                  else
                  {
                     real_t temp = Bj[i];
                     Bj[i] = Bj[i+1];
                     Bj[i+1] = temp - DL[i]*Bj[i];
                  }
               }
               Bj[n-1] = Bj[n-1]/D[n-1];
               if (n > 1)
               {
                  Bj[n-2] = (Bj[n-2]-DU[n-2]*Bj[n-1])/D[n-2];
               }
               for (int_t i = n-3; i >= 0; --i)
               {
                  Bj[i] = (Bj[i] - DU[i]*Bj[i+1]-DU2[i]*Bj[i+2])/D[i];
               }
               Bj += ldB;
            }
         }
      }
      else
      {
         if (nrhs == 1)
         {
            B[0] = B[0]/D[0];
            if (n > 1)
            {
               B[1] = (B[1]-DU[0]*B[0])/D[1];
            }
            for (int_t i = 2; i < n; ++i)
            {
               B[i] = (B[i]-DU[i-1]*B[i-1]-DU2[i-2]*B[i-2])/D[i];
            }
            
            for (int_t i = n-2; i >= 0; --i)
            {
               int_t ip = ipiv[i];
               real_t temp = B[i]-DL[i]*B[i+1];
               B[i] = B[ip];
               B[ip] = temp;
            }
         }
         else
         {
            real_t * Bj = B;
            for (int_t j = 0; j < nrhs; ++j)
            {
               Bj[0] = Bj[0]/D[0];
               if ( n > 1)
                  Bj[1] = (Bj[1]-DU[0]*Bj[0])/D[1];
               for (int_t i = 2; i < n; ++i)
               {
                  Bj[i] = (Bj[i]-DU[i-1]*Bj[i-1]-DU2[i-2]*Bj[i-2])/D[i];
               }
               for (int_t i = n-2; i >= 0; --i)
               {
                  if (ipiv[i] == i)
                  {
                     Bj[i] -= DL[i]*Bj[i+1];
                  }
                  else
                  {
                     real_t temp = Bj[i+1];
                     Bj[i+1] = Bj[i]-DL[i]*temp;
                     Bj[i] = temp;
                  }
               }
               Bj += ldB;
            }
            
         }
      }
      return 0;
   }
   
   /// @brief Solves a complex system of equations of the form
   ///
   ///      A * X = B      or    A' * X = B     or    A.' * X = B
   ///
   /// using the LU factorization of the tridiagonal matrix A computed by GTTRF.
   /// @return 0 if success
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param itrans Specifies the form of the system of equations.
   ///
   ///      0 = no transpose
   ///      1 = transpose
   ///      2 = conjugate transpose
   ///
   /// @param n Order of the matrix A.  n >= 0
   /// @param nrhs Number of columns of the matrix B.  nrhs >= 0
   /// @param DL Complex array, size n-1.  On entry, should contain the (n-1) multipliers that define the matrix L from the LU factorization of A.
   /// @param D Complex array, size n.  On entry, should contain the n diagonal elements of the upper triangular matrix U from the LU factorization of A.
   /// @param DU Complex array, size n-1.  On entry, the (n-1) elements of the first superdiagonal of U.
   /// @param DU2 Complex array, size n-2.  On entry, the (n-2) elements of the second superdiagonal of U.
   /// @param ipiv Integer array, size n.  The pivot indices from the factorization LU.
   /// @param B Complex matrix size ldB-by-nrhs.  On exit, the solution matrix X.
   /// @param ldB Column length of the matrix B.
   /// @ingroup COMP
   
   template< typename real_t>
   int_t GTTS2(const int_t itrans, const int_t n, const int_t nrhs, complex<real_t> * const DL, complex<real_t> * const D, complex<real_t> * const DU, complex<real_t> * const DU2, int_t * const ipiv, complex<real_t> * const B, const int_t ldB)
   {
      if (itrans != 0 && itrans != 1 && itrans != 2)
         return -1;
      if (n < 0)
         return -2;
      if (nrhs < 0)
         return -3;
      if (ldB < n)
         return -10;
      
      if (n == 0 || nrhs == 0)
         return 0;
      
      complex<real_t> temp;
      if (itrans == 0)
      {
         if (nrhs == 1)
         {
            for (int_t i = 0; i < n-1; ++i)
            {
               if (ipiv[i] == i)
               {
                  B[i+1] -= DL[i]*B[i];
               }
               else
               {
                  temp = B[i];
                  B[i] = B[i+1];
                  B[i+1] = temp - DL[i]*B[i];
               }
            }
            B[n-1] = B[n-1]/D[n-1];
            if (n > 1)
            {
               B[n-2] = (B[n-2]-DU[n-2]*B[n-1])/D[n-2];
            }
            for (int_t i = n-3; i >= 0; --i)
            {
               B[i] = (B[i]-DU[i]*B[i+1]-DU2[i]*B[i+2])/D[i];
            }
         }
         else
         {
            complex<real_t> * Bj = B;
            for (int_t j = 0; j < nrhs; ++j)
            {
               for (int_t i = 0; i < n-1; ++i)
               {
                  if (ipiv[i] == i)
                     Bj[i+1] -= DL[i]*Bj[i];
                  else
                  {
                     temp = Bj[i];
                     Bj[i] = Bj[i+1];
                     Bj[i+1] = temp - DL[i]*Bj[i];
                  }
               }
               Bj[n-1] = Bj[n-1]/D[n-1];
               if (n > 1)
               {
                  Bj[n-2] = (Bj[n-2]-DU[n-2]*Bj[n-1])/D[n-2];
               }
               for (int_t i = n-3; i >= 0; --i)
               {
                  Bj[i] = (Bj[i] - DU[i]*Bj[i+1]-DU2[i]*Bj[i+2])/D[i];
               }
               Bj += ldB;
            }
         }
      }
      else if (itrans == 1)
      {
         if (nrhs == 1)
         {
            B[0] = B[0]/D[0];
            if (n > 1)
            {
               B[1] = (B[1]-DU[0]*B[0])/D[1];
            }
            for (int_t i = 2; i < n; ++i)
            {
               B[i] = (B[i]-DU[i-1]*B[i-1]-DU2[i-2]*B[i-2])/D[i];
            }
            
            for (int_t i = n-2; i >= 0; --i)
            {
               if (ipiv[i] == i)
               {
                  B[i] -= DL[i]*B[i+1];
               }
               else
               {
                  temp = B[i+1];
                  B[i+1] = B[i] - DL[i]*temp;
                  B[i] = temp;
               }
            }
         }
         else
         {
            complex<real_t> * Bj = B;
            for (int_t j = 0; j < nrhs; ++j)
            {
               Bj[0] = Bj[0]/D[0];
               if ( n > 1)
                  Bj[1] = (Bj[1]-DU[0]*Bj[0])/D[1];
               for (int_t i = 2; i < n; ++i)
               {
                  Bj[i] = (Bj[i]-DU[i-1]*Bj[i-1]-DU2[i-2]*Bj[i-2])/D[i];
               }
               for (int_t i = n-2; i >= 0; --i)
               {
                  if (ipiv[i] == i)
                  {
                     Bj[i] -= DL[i]*Bj[i+1];
                  }
                  else
                  {
                     temp = Bj[i+1];
                     Bj[i+1] = Bj[i]-DL[i]*temp;
                     Bj[i] = temp;
                  }
               }
               Bj += ldB;
            }
         }
      }
      else 
      {
         if (nrhs == 1)
         {
            B[0] = B[0]/conj(D[0]);
            if (n > 1)
            {
               B[1] = (B[1]-conj(DU[0])*B[0])/conj(D[1]);
            }
            for (int_t i = 2; i < n; ++i)
            {
               B[i] = (B[i]-conj(DU[i-1])*B[i-1]-conj(DU2[i-2])*B[i-2])/conj(D[i]);
            }
            
            for (int_t i = n-2; i >= 0; --i)
            {
               if (ipiv[i] == i)
               {
                  B[i] -= conj(DL[i])*B[i+1];
               }
               else
               {
                  temp = B[i+1];
                  B[i+1] = B[i] - conj(DL[i])*temp;
                  B[i] = temp;
               }
            }
         }
         else
         {
            complex<real_t> * Bj = B;
            for (int_t j = 0; j < nrhs; ++j)
            {
               Bj[0] = Bj[0]/conj(D[0]);
               if ( n > 1)
                  Bj[1] = (Bj[1]-conj(DU[0])*Bj[0])/conj(D[1]);
               for (int_t i = 2; i < n; ++i)
               {
                  Bj[i] = (Bj[i]-conj(DU[i-1])*Bj[i-1]-conj(DU2[i-2])*Bj[i-2])/conj(D[i]);
               }
               for (int_t i = n-2; i >= 0; --i)
               {
                  if (ipiv[i] == i)
                  {
                     Bj[i] -= conj(DL[i])*Bj[i+1];
                  }
                  else
                  {
                     temp = Bj[i+1];
                     Bj[i+1] = Bj[i]-conj(DL[i])*temp;
                     Bj[i] = temp;
                  }
               }
               Bj += ldB;
            }
         }
      }
      return 0;
   }

}

#endif
