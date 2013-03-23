//
//  print.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/23/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _print_h
#define _print_h

/// @file print.h Prints a matrix to standard output.

#include <cctype>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <latl.h>

namespace LATL
{
   /// @brief Prints a real rectangular matrix to standard output.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if ith parameter is invalid.
   /// @param m Number of rows in matrix.
   /// @param n Number of columns in matrix.
   /// @param A Real matrix of size m-by-n.
   /// @param ldA Column length of matrix A.
   /// @param[in] outfile Output stream to write matrix.  Defaults to standard output.

   template <typename real_t>
   int PRINT(int m, int n, real_t *A, int ldA,std::ostream &outfile=std::cout)
   {
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;
      using std::endl;
      using std::setprecision;
      using std::setw;
      using std::numeric_limits;
      using std::floor;
#ifdef MPREAL
      int prec=floor((numeric_limits<real_t>::digits())/3.32193);
#else
      int prec=floor((numeric_limits<real_t>::digits)/3.32193);
#endif
      outfile << setprecision(prec);

      for(int i=0;i<m;i++)
      {
         for(int j=0;j<n;j++)
            outfile << setw(prec+7) << A[i+j*ldA] << "\t";
         outfile << endl;
      }
      outfile << setprecision(6);
      
      return 0;
   }

   /// @brief Prints a complex rectangular matrix to standard output.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if ith parameter is invalid.
   /// @param m Number of rows in matrix.
   /// @param n Number of columns in matrix.
   /// @param A Complex matrix of size m-by-n.
   /// @param ldA Column length of matrix A.
   /// @param[in] outfile Output stream to write matrix.  Defaults to standard output.

   template <typename real_t>
   int PRINT(int m, int n, complex<real_t> *A, int ldA,std::ostream &outfile=std::cout)
   {
      if(m<0)
         return -1;
      else if(n<0)
         return -2;
      else if(ldA<m)
         return -4;
      using std::endl;
      using std::setprecision;
      using std::setw;
      using std::numeric_limits;
      using std::floor;

#ifdef MPREAL
      int prec=floor((numeric_limits<real_t>::digits())/3.32193);
#else
      int prec=floor((numeric_limits<real_t>::digits)/3.32193);
#endif
      outfile << setprecision(prec);

      for(int i=0;i<m;i++)
      {
         for(int j=0;j<n;j++)
            outfile <<  setw(prec+7) << A[i+j*ldA] << "\t";
         outfile << endl;
      }
      outfile << setprecision(6);

      return 0;
   }

   /// @brief Prints a real triangular matrix to standard output.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if ith parameter is invalid.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param diag Specifies whether or not A is unit triangular as follows:
   ///
   ///        if diag = 'U' or 'u' then A is assumed to be unit triangular
   ///        if diag = 'N' or 'n' then A is not assumed to be unit triangular.
   /// @param n Order of triangular matrix.
   /// @param A Real matrix of size m-by-n.
   /// @param ldA Column length of matrix A.
   /// @param[in] outfile Output stream to write matrix.  Defaults to standard output.

   template <typename real_t>
   int PRINT(char uplo, char diag, int n, real_t *A, int ldA,std::ostream &outfile=std::cout)
   {
      const real_t zero=0.0;
      const real_t one=1.0;
      using std::toupper;
      uplo=toupper(uplo);
      diag=toupper(diag);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -3;
      else if(ldA<n)
         return -5;
      using std::endl;
      using std::setprecision;
      using std::setw;
      using std::numeric_limits;
      using std::floor;

#ifdef MPREAL
      int prec=floor((numeric_limits<real_t>::digits())/3.32193);
#else
      int prec=floor((numeric_limits<real_t>::digits)/3.32193);
#endif
      outfile << setprecision(prec);

      if(uplo=='U')
      {
         for(int i=0;i<n;i++)
         {
            for(int j=0;j<i;j++)
               outfile << setw(prec+7) << zero << "\t";
            if(diag=='U')
               outfile << setw(prec+7) << one << "\t";
            else
               outfile << setw(prec+7) << A[i+i*ldA] << "\t";
            for(int j=i+1;j<n;j++)
               outfile << setw(prec+7) << A[i+j*ldA] << "\t";
            outfile << endl;
         }
      }
      else
      {
         for(int i=0;i<n;i++)
         {
            for(int j=0;j<i;j++)
               outfile << setw(prec+7) << A[i+j*ldA] << "\t";
            if(diag=='U')
               outfile << setw(prec+7) << one << "\t";
            else
               outfile << setw(prec+7) << A[i+i*ldA] << "\t";
            for(int j=i+1;j<n;j++)
               outfile << setw(prec+7) << zero << "\t";
            outfile << endl;
         }
      }
      
      outfile << setprecision(6);
      
      return 0;
   }

   /// @brief Prints a complex triangular matrix to standard output.
   /// @tparam real_t Floating point type.
   /// @return 0 if success.
   /// @return -i if ith parameter is invalid.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param diag Specifies whether or not A is unit triangular as follows:
   ///
   ///        if diag = 'U' or 'u' then A is assumed to be unit triangular
   ///        if diag = 'N' or 'n' then A is not assumed to be unit triangular.
   /// @param n Order of triangular matrix.
   /// @param A Complex matrix of size m-by-n.
   /// @param ldA Column length of matrix A.
   /// @param[in] outfile Output stream to write matrix.  Defaults to standard output.

   template <typename real_t>
   int PRINT(char uplo, char diag, int n, complex<real_t> *A, int ldA,std::ostream &outfile=std::cout)
   {
      const real_t zero=0.0;
      const real_t one=1.0;
      using std::toupper;
      uplo=toupper(uplo);
      diag=toupper(diag);
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if(n<0)
         return -3;
      else if(ldA<n)
         return -5;
      using std::endl;
      using std::setprecision;
      using std::setw;
      using std::numeric_limits;
      using std::floor;

#ifdef MPREAL
      int prec=floor((numeric_limits<real_t>::digits())/3.32193);
#else
      int prec=floor((numeric_limits<real_t>::digits)/3.32193);
#endif
      outfile << setprecision(prec);

      if(uplo=='U')
      {
         for(int i=0;i<n;i++)
         {
            for(int j=0;j<i;j++)
               outfile << setw(prec+7) << zero << "\t";
            if(diag=='U')
               outfile << setw(prec+7) << one << "\t";
            else
               outfile << setw(prec+7) << A[i+i*ldA] << "\t";
            for(int j=i+1;j<n;j++)
               outfile << setw(prec+7) << A[i+j*ldA] << "\t";
            outfile << endl;
         }
      }
      else
      {
         for(int i=0;i<n;i++)
         {
            for(int j=0;j<i;j++)
               outfile << setw(prec+7) << A[i+j*ldA] << "\t";
            if(diag=='U')
               outfile << setw(prec+7) << one << "\t";
            else
               outfile << setw(prec+7) << A[i+i*ldA] << "\t";
            for(int j=i+1;j<n;j++)
               outfile << setw(prec+7) << zero << "\t";
            outfile << endl;
         }
      }

      outfile << setprecision(6);

      return 0;
   }
}

#endif
