//
//  dlamch.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 1/1/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "lapack.h"
#include <limits>
#include <cctype>

using std::numeric_limits;
using std::toupper;

double dlamch_(char &OPT)
{
   const double one=1.0;
   const double zero=0.0;
   const double prec=numeric_limits<double>::epsilon();
   const double rnd=numeric_limits<double>::round_error();
   const double eps=prec*rnd;
   const double maxval=numeric_limits<double>::max();
   const double minval=numeric_limits<double>::min();
   const double small=one/maxval;
   const double sfmin=(small>minval) ? small*(one+eps) : minval;
   const double base=numeric_limits<double>::radix;
   const double num=numeric_limits<double>::digits;
   const double round=numeric_limits<double>::round_style;
   const double min_exp=numeric_limits<double>::min_exponent;
   const double underflow=numeric_limits<double>::min();
   const double max_exp=numeric_limits<double>::max_exponent;
   const double overflow=numeric_limits<double>::max();
   
   const char opt=toupper(OPT);
   
   if(opt=='E')
      return eps;
   else if(opt=='S')
      return sfmin;
   else if(opt=='B')
      return base;
   else if(opt=='P')
      return prec;
   else if(opt=='N')
      return num;
   else if(opt=='R')
      return round;
   else if(opt=='M')
      return min_exp;
   else if(opt=='U')
      return underflow;
   else if(opt=='L')
      return max_exp;
   else if(opt=='O')
      return overflow;
   else
      return zero;
}


