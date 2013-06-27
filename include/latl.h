//
//  latl.h
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/15/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#ifndef _latl_h
#define _latl_h

/// @file latl.h Global types used in LATL.

/// @mainpage
/// @section intro Introduction
/// The Linear Algebra Template Library is a precision neutral dense linear algebra
/// library which implements the core functionality of
/// <A HREF="http://netlib.org/lapack/">LAPACK</A>
/// in C++ using function templates.  The floating point type is
/// specified by a template parameter, enabling arbitrary precision computations
/// by using a multiprecision class such as Pavel Holoborodko's
/// <A HREF="http://www.holoborodko.com/pavel/mpfr/">mpfr::mpreal</A> or
/// Christian Schneider's
/// <A HREF="http://chschneider.eu/programming/mpfr_real/">mpfr::real</A>.
/// Matrices are implemented as pointers to column-major contiguous arrays,
/// rather than using objects, in order to simplify operability with other
/// libraries as well as preserve compatibility with LAPACK.
///
/// @section code Download
/// Source code and git access is available on GitHub at
/// <A HREF="http://github.com/langou/latl">http://github.com/langou/latl</A>.
///
/// @defgroup DRIV Driver Routines
/// @defgroup COMP Computational Routines
/// @defgroup AUX Auxiliary Routines
/// @defgroup BLAS BLAS Routines
/// @defgroup MATGEN Matrix Generation Routines

#include <complex>
using std::complex;
#include <cstddef>

namespace LATL
{
   /// Integer type used in LATL namespace.
   
   typedef std::ptrdiff_t int_t;
}

#endif
