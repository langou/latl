//
//  tpmv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/13/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _tpmv_h
#define _tpmv_h

/// @file tpmv.h Performs triangular matrix-vector multiplication using packed storage.

#include <cctype>
#include "latl.h"

namespace LATL
{
   
   /// @brief Performs real triangular matrix-vector multiplication using packed storage.
   /// 
   /// For a real upper or lower triangular matrix A and real vector x,
   ///
   ///        x := A * x  or  x := A' * x
   ///
   /// is computed, where A uses packed storage format. 
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param trans Specifies wheather the transpose of A is to be used or not:
   ///
   ///        if trans = 'N' or 'n' then x := A*x
   ///        if trans = 'T' or 't' then x := A'*x
   ///        if trans = 'C' or 'c' then x := A'*x
   ///
   /// @param diag specifies whether or not A is unit triangular as follows:
   ///
   ///        if diag = 'U' or 'u' then A is assumed to be unit triangular
   ///        if diag = 'N' or 'n' then A is not assumed to be unit triangular.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param A Pointer to packed real triangular n-by-n matrix A.  The packed storage format stores a triangular matrix
   /// using variable length columns so that no space is wasted.  If A is upper triangular, the first column consists only
   /// of the diagonal element, which has length 1.  The second column has length 2, the third length 3, and so on until the
   /// nth column, which has length n.  Similarly, if A is lower triangular, the first column has length n, the second length n-1,
   /// and so on until the nth column, which has length 1.  The entire n-by-n triangular matrix is stored using n(n+1)/2 elements.
   /// @param x Pointer to real vector x.
   /// @param incx Increment of the vector x.  incx!=0
   /// @ingroup MATV

   template <typename real_t>
   int TPMV(char uplo, char trans, char diag, int_t n, real_t *A, real_t *x, int_t incx)
   {
      using std::toupper;

      int_t i,j,kx,jx,ix;

      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -2;
      else if((diag!='U')&&(diag!='N'))
         return -3;
      else if(n<0)
         return -4;
      else if(incx==0)
         return -7;
      else if(n==0)
         return 0;
      
      bool nounit=(diag=='N');
      
      if(incx==1)
      {
         if(trans=='N')
         {
            if(uplo=='U')
            {  
               for(j=0;j<n;j++)
               {
                  for(i=0;i<j;i++)
                     x[i]+=x[j]*A[i];
                  if(nounit)
                     x[j]*=A[j];
                  A+=j+1;
               }
            }
            else
            {
               A+=n*(n+1)/2;
               for(j=n-1;j>=0;j--)
               {
                  A-=n-j;
                  for(i=n-1;i>j;i--)
                     x[i]+=x[j]*A[i-j];
                  if(nounit)
                     x[j]*=A[0];
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               A+=n*(n+1)/2;
               for(j=n-1;j>=0;j--)
               {
                  A-=j+1;
                  if(nounit)
                     x[j]*=A[j];
                  for(i=j-1;i>=0;i--)
                     x[j]+=A[i]*x[i];
               }
            }
            else
            {
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[j]*=A[0];
                  for(i=j+1;i<n;i++)
                     x[j]+=A[i-j]*x[i];
                  A+=n-j;
               }
            }
         }
      }
      else
      {
         kx=(incx>0)?0:(1-n)*incx;
         if(trans=='N')
         {
            if(uplo=='U')
            {  
               jx=kx;
               for(j=0;j<n;j++)
               {
                  ix=kx;
                  for(i=0;i<j;i++)
                  {
                     x[ix]+=x[jx]*A[i];
                     ix+=incx;
                  }
                  if(nounit)
                     x[jx]*=A[j];
                  A+=j+1;
                  jx+=incx;
               }
            }
            else
            {
               A+=n*(n+1)/2;
               kx+=n*incx;
               jx=kx;
               for(j=n-1;j>=0;j--)
               {
                  A-=n-j;
                  jx-=incx;
                  ix=kx;
                  for(i=n-1;i>j;i--)
                  {
                     ix-=incx;
                     x[ix]+=x[jx]*A[i-j];
                  }
                  if(nounit)
                     x[jx]*=A[0];
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               A+=n*(n+1)/2;
               jx=kx+n*incx;
               for(j=n-1;j>=0;j--)
               {
                  A-=j+1;
                  jx-=incx;
                  ix=jx;
                  if(nounit)
                     x[jx]*=A[j];
                  for(i=j-1;i>=0;i--)
                  {
                     ix-=incx;
                     x[jx]+=A[i]*x[ix];
                  }
               }
            }
            else
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  ix=jx;
                  if(nounit)
                     x[jx]*=A[0];
                  for(i=j+1;i<n;i++)
                  {
                     ix+=incx;
                     x[jx]+=A[i-j]*x[ix];
                  }
                  A+=n-j;
                  jx+=incx;
               }
            }
         }      
      }
      
      return 0;
   }

   /// @brief Performs complex triangular matrix-vector multiplication using packed storage.
   /// 
   /// For a real upper or lower triangular matrix A and real vector x,
   ///
   ///        x := A * x,  x := A.' * x  or  x := A' * x
   ///
   /// is computed, where A uses packed storage format. 
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param trans Specifies wheather the transpose, conjugate transpose or neither is used:
   ///
   ///        if trans = 'N' or 'n' then x := A*x
   ///        if trans = 'T' or 't' then x := A.'*x
   ///        if trans = 'C' or 'c' then x := A'*x
   ///
   /// @param diag specifies whether or not A is unit triangular as follows:
   ///
   ///        if diag = 'U' or 'u' then A is assumed to be unit triangular
   ///        if diag = 'N' or 'n' then A is not assumed to be unit triangular.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param A Pointer to packed complex triangular n-by-n matrix A.  The packed storage format stores a triangular matrix
   /// using variable length columns so that no space is wasted.  If A is upper triangular, the first column consists only
   /// of the diagonal element, which has length 1.  The second column has length 2, the third length 3, and so on until the
   /// nth column, which has length n.  Similarly, if A is lower triangular, the first column has length n, the second length n-1,
   /// and so on until the nth column, which has length 1.  The entire n-by-n triangular matrix is stored using n(n+1)/2 elements.
   /// @param x Pointer to complex vector x.
   /// @param incx Increment of the vector x.  incx!=0
   /// @ingroup MATV

   template <typename real_t>
   int TPMV(char uplo, char trans, char diag, int_t n, complex<real_t> *A, complex<real_t> *x, int_t incx)
   {
      using std::conj;
      using std::toupper;
      
      int_t i,j,kx,jx,ix;

      uplo=toupper(uplo);
      trans=toupper(trans);
      diag=toupper(diag);
      
      if((uplo!='U')&&(uplo!='L'))
         return -1;
      else if((trans!='N')&&(trans!='T')&&(trans!='C'))
         return -2;
      else if((diag!='U')&&(diag!='N'))
         return -3;
      else if(n<0)
         return -4;
      else if(incx==0)
         return -7;
      else if(n==0)
         return 0;
      
      bool nounit=(diag=='N');
      
      if(incx==1)
      {
         if(trans=='N')
         {
            if(uplo=='U')
            {  
               for(j=0;j<n;j++)
               {
                  for(i=0;i<j;i++)
                     x[i]+=x[j]*A[i];
                  if(nounit)
                     x[j]*=A[j];
                  A+=j+1;
               }
            }
            else
            {
               A+=n*(n+1)/2;
               for(j=n-1;j>=0;j--)
               {
                  A-=n-j;
                  for(i=n-1;i>j;i--)
                     x[i]+=x[j]*A[i-j];
                  if(nounit)
                     x[j]*=A[0];
               }
            }
         }
         else if(trans=='T')
         {
            if(uplo=='U')
            {
               A+=n*(n+1)/2;
               for(j=n-1;j>=0;j--)
               {
                  A-=j+1;
                  if(nounit)
                     x[j]*=A[j];
                  for(i=j-1;i>=0;i--)
                     x[j]+=A[i]*x[i];
               }
            }
            else
            {
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[j]*=A[0];
                  for(i=j+1;i<n;i++)
                     x[j]+=A[i-j]*x[i];
                  A+=n-j;
               }
            }
         }
         else 
         {
            if(uplo=='U')
            {
               A+=n*(n+1)/2;
               for(j=n-1;j>=0;j--)
               {
                  A-=j+1;
                  if(nounit)
                     x[j]*=conj(A[j]);
                  for(i=j-1;i>=0;i--)
                     x[j]+=conj(A[i])*x[i];
               }
            }
            else
            {
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[j]*=conj(A[0]);
                  for(i=j+1;i<n;i++)
                     x[j]+=conj(A[i-j])*x[i];
                  A+=n-j;
               }
            }
         }
         
      }
      else
      {
         kx=(incx>0)?0:(1-n)*incx;
         if(trans=='N')
         {
            if(uplo=='U')
            {  
               jx=kx;
               for(j=0;j<n;j++)
               {
                  ix=kx;
                  for(i=0;i<j;i++)
                  {
                     x[ix]+=x[jx]*A[i];
                     ix+=incx;
                  }
                  if(nounit)
                     x[jx]*=A[j];
                  A+=j+1;
                  jx+=incx;
               }
            }
            else
            {
               A+=n*(n+1)/2;
               kx+=n*incx;
               jx=kx;
               for(j=n-1;j>=0;j--)
               {
                  A-=n-j;
                  jx-=incx;
                  ix=kx;
                  for(i=n-1;i>j;i--)
                  {
                     ix-=incx;
                     x[ix]+=x[jx]*A[i-j];
                  }
                  if(nounit)
                     x[jx]*=A[0];
               }
            }
         }
         else if(trans=='T')
         {
            if(uplo=='U')
            {
               A+=n*(n+1)/2;
               jx=kx+n*incx;
               for(j=n-1;j>=0;j--)
               {
                  A-=j+1;
                  jx-=incx;
                  ix=jx;
                  if(nounit)
                     x[jx]*=A[j];
                  for(i=j-1;i>=0;i--)
                  {
                     ix-=incx;
                     x[jx]+=A[i]*x[ix];
                  }
               }
            }
            else
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  ix=jx;
                  if(nounit)
                     x[jx]*=A[0];
                  for(i=j+1;i<n;i++)
                  {
                     ix+=incx;
                     x[jx]+=A[i-j]*x[ix];
                  }
                  A+=n-j;
                  jx+=incx;
               }
            }
         }      
         else
         {
            if(uplo=='U')
            {
               A+=n*(n+1)/2;
               jx=kx+n*incx;
               for(j=n-1;j>=0;j--)
               {
                  A-=j+1;
                  jx-=incx;
                  ix=jx;
                  if(nounit)
                     x[jx]*=conj(A[j]);
                  for(i=j-1;i>=0;i--)
                  {
                     ix-=incx;
                     x[jx]+=conj(A[i])*x[ix];
                  }
               }
            }
            else
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  ix=jx;
                  if(nounit)
                     x[jx]*=conj(A[0]);
                  for(i=j+1;i<n;i++)
                  {
                     ix+=incx;
                     x[jx]+=conj(A[i-j])*x[ix];
                  }
                  A+=n-j;
                  jx+=incx;
               }
            }
         }      
      }
      
      return 0;
   }
}
#endif
