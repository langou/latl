//
//  trsv.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 12/23/11.
//  Copyright (c) 2011 University of Colorado Denver. All rights reserved.
//

#ifndef _trsv_h
#define _trsv_h

/// @file trsv.h Solves general triangular system of equations.

#include <cctype>
#include "latl.h"

namespace LATL
{
   /// @brief Solves general real triangular system of equations.
   /// 
   /// For a real upper or lower triangular matrix A and real vector x,
   ///
   ///       A*x=b  or  A'*x=b
   ///
   /// is solved for x.  
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
   /// If A is unit triangular, the diagonal elements are assumed to be unity and are not referenced.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param A Pointer to real triangular n-by-n matrix A.  If uplo = 'U' or 'u, A is upper triangular and the lower triangular
   /// part is not referenced.  If uplo = 'L' or 'l A is lower triangular and the upper triangular part is not referenced.
   /// @param ldA Column length of the matrix A.  ldA>=n
   /// @param x Pointer to real vector x.  On entry, x contains the right hand side vector b.
   /// On exit, x is overwritten with the solution vector x.
   /// @param incx Increment of the vector x.  incx!=0
   /// @ingroup SOLV

   template <typename real_t>
   int TRSV(char uplo, char trans, char diag, int_t n, real_t *A, int_t ldA, real_t *x, int_t incx)
   {
      using std::toupper;
      int_t i,j,ix,jx,kx;
      bool nounit;

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
      else if(ldA<n)
         return -6;
      else if(incx==0)
         return -8;
      else if(n==0)
         return 0;
      
      nounit=(diag=='N');
      
      if(incx==1)
      {
         if(trans=='N')
         {
            if(uplo=='U')
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  if(nounit)
                     x[j]=x[j]/A[j];
                  for(i=j-1;i>=0;i--)
                     x[i]-=x[j]*A[i];
               }
            }
            else
            {
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[j]=x[j]/A[j];
                  for(i=j+1;i<n;i++)
                     x[i]-=x[j]*A[i];
                  A+=ldA;
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               for(j=0;j<n;j++)
               {
                  for(i=0;i<j;i++)
                     x[j]-=x[i]*A[i];
                  if(nounit)
                     x[j]=x[j]/A[j];
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  for(i=n-1;i>j;i--)
                     x[j]-=x[i]*A[i];
                  if(nounit)
                     x[j]=x[j]/A[j];
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
               A+=n*ldA;
               jx=kx+n*incx;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  jx-=incx;
                  ix=jx;
                  if(nounit)
                     x[jx]=x[jx]/A[j];
                  for(i=j-1;i>=0;i--)
                  {
                     ix-=incx;
                     x[ix]-=x[jx]*A[i];
                  }
               }
            }
            else
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[jx]=x[jx]/A[j];
                  ix=jx;
                  for(i=j+1;i<n;i++)
                  {
                     ix+=incx;
                     x[ix]-=x[jx]*A[i];
                  }
                  jx+=incx;
                  A+=ldA;
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  ix=kx;
                  for(i=0;i<j;i++)
                  {
                     x[jx]-=x[ix]*A[i];
                     ix+=incx;
                  }
                  if(nounit)
                     x[jx]=x[jx]/A[j];
                  jx+=incx;
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               kx+=n*incx;
               jx=kx;
               for(j=n-1;j>=0;j--)
               {
                  jx-=incx;
                  ix=kx;
                  A-=ldA;
                  for(i=n-1;i>j;i--)
                  {
                     ix-=incx;
                     x[jx]-=x[ix]*A[i];
                  }
                  if(nounit)
                     x[jx]=x[jx]/A[j];
               }
            }
         }
      }
      return 0;
   }
   
   /// @brief Solves general complex triangular system of equations.
   /// 
   /// For a complex upper or lower triangular matrix A and complex vector x,
   ///
   ///       A*x=b  or  A'*x=b  or A.'*x=b
   ///
   /// is solved for x.  
   /// @return 0 if success.
   /// @return -i if the ith argument is invalid.
   /// @tparam real_t Floating point type.
   /// @param uplo Specifies whether A is stored as upper or lower triangular.
   ///
   ///        if uplo = 'U' or 'u' then A is upper triangular
   ///        if uplo = 'L' or 'l' then A is lower triangular
   /// @param trans Specifies wheather the transpose or conjugate transpose of A is to be used:
   ///
   ///        if trans = 'N' or 'n' then x := A*x
   ///        if trans = 'T' or 't' then x := A.'*x
   ///        if trans = 'C' or 'c' then x := A'*x
   ///
   /// @param diag specifies whether or not A is unit triangular as follows:
   ///
   ///        if diag = 'U' or 'u' then A is assumed to be unit triangular
   ///        if diag = 'N' or 'n' then A is not assumed to be unit triangular.
   /// If A is unit triangular, the diagonal elements are assumed to be unity and are not referenced.
   /// @param n Specifies the order of the matrix A.  n>=0
   /// @param A Pointer to complex triangular n-by-n matrix A.  If uplo = 'U' or 'u, A is upper triangular and the lower triangular
   /// part is not referenced.  If uplo = 'L' or 'l A is lower triangular and the upper triangular part is not referenced.
   /// @param ldA Column length of the matrix A.  ldA>=n
   /// @param x Pointer to complex vector x.  On entry, x contains the right hand side vector b.
   /// On exit, x is overwritten with the solution vector x.
   /// @param incx Increment of the vector x.  incx!=0
   /// @ingroup SOLV

   template <typename real_t>
   int TRSV(char uplo, char trans, char diag, int_t n, complex<real_t> *A, int_t ldA, complex<real_t> *x, int_t incx)
   {
      using std::conj;
      using std::toupper;
      int_t i,j,ix,jx,kx;
      bool nounit;
      
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
      else if(ldA<n)
         return -6;
      else if(incx==0)
         return -8;
      else if(n==0)
         return 0;
      
      nounit=(diag=='N');
      
      if(incx==1)
      {
         if(trans=='N')
         {
            if(uplo=='U')
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  if(nounit)
                     x[j]=x[j]/A[j];
                  for(i=j-1;i>=0;i--)
                     x[i]-=x[j]*A[i];
               }
            }
            else
            {
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[j]=x[j]/A[j];
                  for(i=j+1;i<n;i++)
                     x[i]-=x[j]*A[i];
                  A+=ldA;
               }
            }
         }
         else if(trans=='T')
         {
            if(uplo=='U')
            {
               for(j=0;j<n;j++)
               {
                  for(i=0;i<j;i++)
                     x[j]-=x[i]*A[i];
                  if(nounit)
                     x[j]=x[j]/A[j];
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  for(i=n-1;i>j;i--)
                     x[j]-=x[i]*A[i];
                  if(nounit)
                     x[j]=x[j]/A[j];
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               for(j=0;j<n;j++)
               {
                  for(i=0;i<j;i++)
                     x[j]-=x[i]*conj(A[i]);
                  if(nounit)
                     x[j]=x[j]/conj(A[j]);
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  for(i=n-1;i>j;i--)
                     x[j]-=x[i]*conj(A[i]);
                  if(nounit)
                     x[j]=x[j]/conj(A[j]);
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
               A+=n*ldA;
               jx=kx+n*incx;
               for(j=n-1;j>=0;j--)
               {
                  A-=ldA;
                  jx-=incx;
                  ix=jx;
                  if(nounit)
                     x[jx]=x[jx]/A[j];
                  for(i=j-1;i>=0;i--)
                  {
                     ix-=incx;
                     x[ix]-=x[jx]*A[i];
                  }
               }
            }
            else
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  if(nounit)
                     x[jx]=x[jx]/A[j];
                  ix=jx;
                  for(i=j+1;i<n;i++)
                  {
                     ix+=incx;
                     x[ix]-=x[jx]*A[i];
                  }
                  jx+=incx;
                  A+=ldA;
               }
            }
         }
         else if(trans=='T')
         {
            if(uplo=='U')
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  ix=kx;
                  for(i=0;i<j;i++)
                  {
                     x[jx]-=x[ix]*A[i];
                     ix+=incx;
                  }
                  if(nounit)
                     x[jx]=x[jx]/A[j];
                  jx+=incx;
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               kx+=n*incx;
               jx=kx;
               for(j=n-1;j>=0;j--)
               {
                  jx-=incx;
                  ix=kx;
                  A-=ldA;
                  for(i=n-1;i>j;i--)
                  {
                     ix-=incx;
                     x[jx]-=x[ix]*A[i];
                  }
                  if(nounit)
                     x[jx]=x[jx]/A[j];
               }
            }
         }
         else
         {
            if(uplo=='U')
            {
               jx=kx;
               for(j=0;j<n;j++)
               {
                  ix=kx;
                  for(i=0;i<j;i++)
                  {
                     x[jx]-=x[ix]*conj(A[i]);
                     ix+=incx;
                  }
                  if(nounit)
                     x[jx]=x[jx]/conj(A[j]);
                  jx+=incx;
                  A+=ldA;
               }
            }
            else
            {
               A+=n*ldA;
               kx+=n*incx;
               jx=kx;
               for(j=n-1;j>=0;j--)
               {
                  jx-=incx;
                  ix=kx;
                  A-=ldA;
                  for(i=n-1;i>j;i--)
                  {
                     ix-=incx;
                     x[jx]-=x[ix]*conj(A[i]);
                  }
                  if(nounit)
                     x[jx]=x[jx]/conj(A[j]);
               }
            }
         }      
      }
      return 0;
   }
}
#endif

