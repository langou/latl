//
//  timer.cpp
//  Linear Algebra Template Library
//
//  Created by Rodney James on 3/25/13.
//  Copyright (c) 2013 University of Colorado Denver. All rights reserved.
//

#include "timer.h"

#if defined(__hpux)
#include <time.h>
#include <stdio.h>

double LATL::Timer::Clock()
{
  double t;
  struct timeval buffer;
  struct timezone dummy;
  gettimeofday (&buffer, &dummy);
  t = (double)(buffer.tv_sec*1000000 + buffer.tv_usec);
  return (t);
}

#elif defined (__sun)
#include <time.h>
double LATL::Timer::Clock()
{
  struct timespec tp;
  clock_gettime(CLOCK_HIGHRES,&tp);
  double usec=1.0e-03*(double)tp.tv_nsec+1.0e+06*(double)tp.tv_sec;
  return(usec);
}

#elif defined(linux) 
#include <sys/time.h>

double LATL::Timer::Clock()
{
  double t;
  struct timeval buffer;
  struct timezone dummy;
  gettimeofday (&buffer, &dummy);
  t = (double)(buffer.tv_sec*1000000 + buffer.tv_usec);
  return (t);
}

#elif defined(__alpha) 
#include <sys/time.h>
double LATL::Timer::Clock()
{
  double t;
  struct timeval buffer;
  struct timezone dummy;
  gettimeofday (&buffer, &dummy);
  t = (double)(buffer.tv_sec*1000000+buffer.tv_usec);
  return (t);
}

#elif defined(__sgi) 
#include <time.h>
double LATL::Timer::Clock()
{
  struct timespec tp;
  clock_gettime(CLOCK_REALTIME,&tp);
  double usec=1.0e-03*(double)tp.tv_nsec+1.0e+06*(double)tp.tv_sec;
  return(usec);
}

#elif defined(_AIX) 
#include <sys/time.h>
#include <sys/systemcfg.h>
double LATL::Timer::Clock()
{
  double t;
  timebasestruct_t T;
  read_real_time(&T,TIMEBASE_SZ);
  time_base_to_time(&T,TIMEBASE_SZ);
  t=(double)(T.tb_low)*1.0e-03+(double)(T.tb_high)*1.0e+6;
  return (t);
}
#elif defined(__APPLE__) 
#include <sys/time.h>
double LATL::Timer::Clock()
{
  double t;
  struct timeval buffer;
  struct timezone dummy;
  gettimeofday (&buffer, &dummy);
  t = (double)(buffer.tv_sec*1000000 + buffer.tv_usec);
  return (t);
}
#else
#include <time.h>
double LATL::Timer::Clock(void)
{
  clock_t t0;
  double c=1.0e+06/(double)CLOCKS_PER_SEC;
  double t;
  t=c*(double)clock();
  return(t);
}
#endif

LATL::Timer::Timer()
{
  Clear();
}

LATL::Timer::~Timer()
{
}

void LATL::Timer::Start()
{
  StopWatch=-Clock();
}

void LATL::Timer::Stop()
{
  StopWatch+=Clock();
  Total+=StopWatch;
}

void LATL::Timer::Clear()
{
  Total=0.0;
  StopWatch=0.0;
}

double LATL::Timer::Time()
{
  return(Total);
}
