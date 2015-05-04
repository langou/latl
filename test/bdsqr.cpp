//
//  bdsqr.cpp
//  Linear Algebra Template Library
//
//  Created by Philippe Theveny on 4/10/15.
//  Copyright (c) 2014 University of Colorado Denver. All rights reserved.
//
//  Check singular values of a bidiagonal matrix and report the relative error
//
// ||B - U * S * V'|| / ||B||
// along with the distances to orthogonality of the computed U and V
//
// ||U*U' - I||, ||U'*U - I||, ||V*V' - I||, and ||V'*V - I||
//

#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <time.h>

#ifdef REAL
# include "real.hpp"
#endif

#ifdef MPREAL
# include "mpreal.h"
#endif

#include <load.h>
#include <print.h>
#include "timer.h"

#include <latl.h>
#include <bdsqr.h>
#include <gemm.h>
#include <lacpy.h>
#include <lamch.h>
#include <lange.h>
#include <laset.h>

namespace TestBDSQR
{
  using namespace LATL;

  class TestParameters;

  class Test {
  public:
    // test parameters
    int_t threshold;
    std::string inLeft;
    std::string outLeft;
    int_t nrLeft;
    int_t ncLeft;
    std::string inBidiag;
    std::string outSingular;
    char uplo;
    int_t n;
    std::string inRight;
    std::string outRight;
    int_t nrRight;
    int_t ncRight;
    std::string inAux;
    std::string outAux;
    int_t nrAux;
    int_t ncAux;
    std::string outLog;

    // statistics
    std::string datestart;
    std::string datestop;
    double time;
    double residual;
    double orthogonality[4];
    int_t info;

    virtual int_t Run () = 0;

    virtual ~Test () {};
  };

  template <typename real_t, typename matrix_t>
  class TestTemplate : public Test {
  private:
    // data
    matrix_t* D;
    matrix_t* E;
    matrix_t* Bidiag;
    matrix_t* Singular;
    matrix_t* Left;
    matrix_t* Right;
    matrix_t* Auxiliary;

    real_t eps;

    // possibly expensive constants
    constexpr static const real_t zero = 0.0;
    constexpr static const real_t one  = 1.0;
    constexpr static const real_t negone = -1.0;

    void
    ReadParameters (std::istream &file)
    {
      using std::cout;
      using std::endl;

      file >>threshold;
      file >>inLeft;
      file >>outLeft;
      file >>nrLeft;
      file >>ncLeft;
      file >>inBidiag;
      file >>outSingular;
      file >>uplo;
      file >>n;
      file >>inRight;
      file >>outRight;
      file >>nrRight;
      file >>ncRight;
      file >>inAux;
      file >>outAux;
      file >>nrAux;
      file >>ncAux;
      file >>outLog;
      if (file.fail())
      {
        cout <<"Error: Cannot read input parameters." <<endl;
        exit (-1);
      }
    }

    void
    ReadBidiagonalMatrix ()
    {
      using std::cout;
      using std::endl;
      using std::ifstream;

      ifstream iFile;
      const int offset = (uplo=='U')? +1 : -1;

      iFile.open(inBidiag);
      for (int_t i=0; i<n; i++)
        iFile >>D[i] >>E[i];
      if (iFile.fail())
      {
        cout <<"Error: Cannot read file '" <<inBidiag <<"'." <<endl;
        exit (-1);
      }
      iFile.close();

      for (int_t i=0; i<n; i++)
      {
        for (int_t j=0; j<n; j++)
        {
          if (i == j)
            Bidiag[i+j*n] = D[i];
          else if (j == i+offset )
            Bidiag[i+j*n] = E[(offset>0)?i:j];
          else
            Bidiag[i+j*n]= zero;
        }
      }
    }

    bool
    HaveSameComponents (int_t m1, int_t n1, real_t *M1, int_t ldM1,
                        int_t m2, int_t n2, real_t *M2, int_t ldM2)
    {
      using std::min;
      using std::abs;

      real_t tol(threshold*eps);
      int_t minm = min(m1, m2);
      int_t minn = min(n1, n2);

      for (int_t j=0; j<minn; j++)
        for (int_t i=0; i<minm; i++)
          if (abs(M1[i+j*ldM1]-M2[i+j*ldM2]) >= tol)
            return false;

      return true;
    }

    bool
    MatchReference (real_t &tol, std::ifstream &ref, int_t &r, int_t &c,
                    matrix_t *A, int_t ldA)
    {
      using std::cout;
      using std::endl;
      real_t comp;

      int_t nr = r;
      int_t nc = c;

      for (r=0; r<nr; r++)
        for (c=0; c<nc; c++)
        {
          ref >>comp;
          if (abs(comp) > zero)
          {
            if (abs(A[r+c*ldA] - comp) >= tol*abs(comp))
            {
              cout <<A[r+c*ldA] <<"=" <<comp <<"*" <<(A[r+c*ldA] - comp)/abs(comp) <<endl;
              return false;
            }
          }
          else if (abs(A[r+c*ldA] - comp) >= tol)
          {
            cout <<A[r+c*ldA] <<"=" <<comp <<"+" <<A[r+c*ldA] - comp <<endl;
            return false;
          }
        }

      return true;
    }

    bool
    IsOrthogonal (real_t &tol, int_t m, int_t n, matrix_t *A, int_t ldA,
                  real_t &ortho1, real_t &ortho2)
    {
      matrix_t *M = new matrix_t[m*m];
      matrix_t *N = (m == n) ? M : new matrix_t[n*n];

      // ||A*A'-I||
      LASET<real_t>('a', m, m, 0, 1, M, m);
      GEMM<real_t>('n', 't', m, m, n, one, A, ldA, A, ldA, negone, M, m);
      ortho1 = LANGE<real_t>('m', m, m, M, m);

      // ||A'*A-I||
      LASET<real_t>('a', n, n, 0, 1, N, n);
      GEMM<real_t>('t', 'n', n, n, m, one, A, ldA, A, ldA, negone, N, n);
      ortho2 = LANGE<real_t>('m', n, n, N, n);

      delete [] M;
      if (m != n)
        delete [] N;

      return ortho1 < tol && ortho2 < tol;
    }

    void
    CheckDiagonalization ()
    {
      using std::cout;
      using std::endl;

      real_t *fullLeft;
      real_t *fullRight;
      real_t *M1 = new real_t[n*n];
      real_t *M2 = new real_t[n*n];
      real_t normBD = LANGE<real_t>('m', n, n, Bidiag, n);

      if (nrLeft != n || ncRight != n)
      {
        // compute all rows of the left orthogonal matrix and
        // all columns of the right orthogonal matrix

        real_t *newD = new real_t[n];
        real_t *newE = new real_t[n];
        fullLeft = new real_t[n*n];
        fullRight = new real_t[n*n];

        int_t offset = (uplo == 'U') ? n : 1;
        for (int_t i=0; i<n-1; i++)
        {
          D[i]=Bidiag[i+i*n];
          E[i]=Bidiag[i+i*n+offset];
        }
        D[n-1]=Bidiag[n*n-1];
        E[n-1]=zero;
        int_t ret = BDSQR<real_t>(uplo, n, n, n, 0, newD, newE,
                                   fullRight, n, fullLeft, n, NULL, 1);
        if (ret!=0)
        {
          cout <<"Error: BDSQR failed. INFO=" <<ret <<endl;
          exit (1);
        }

        if (!HaveSameComponents(n, 1, D, 1, n, 1, newD, 1))
        {
          cout <<"Error: singular values differ when computed under test"
               <<" conditions and when full left and right transformations"
               <<" are computed." <<endl;
          info |= 8;
        }
        if (!HaveSameComponents(nrLeft, ncLeft, Left, nrLeft,
                                n, n, fullLeft, n))
        {
          cout <<"Error: computed left matrices differ when computed"
               <<" under test conditions and when full left and right"
               <<" transformations are computed." <<endl;
          info |= 16;
        }
        if (!HaveSameComponents(nrRight, nrRight, Right, nrRight,
                                n, n, fullRight, n))
        {
          cout <<"Error: computed right matrices differ when computed"
               <<" under test conditions and when full left and right"
               <<" transformations are computed." <<endl;
          info |= 32;
        }
        delete [] newD;
        delete [] newE;
      }
      else
      {
        fullLeft = Left;
        fullRight = Right;
      }

      LACPY<real_t>('A', n, n, Bidiag, n, M2, n);
      for (int_t j=0; j<n; j++)
        for (int_t i=0; i<n; i++)
          M1[i+j*n]=fullLeft[i+j*n]*D[j];
      GEMM<real_t>('n', 'n', n, n, n, one, M1, n, fullRight, n, negone, M2, n);
      residual = LANGE<real_t>('m', n, n, M2, n);

      const real_t tol(threshold * eps * normBD);
      if (residual >= tol * n)
        info |= 1;

      delete [] M1;
      delete [] M2;
      if (nrLeft != n || ncRight != n)
      {
        delete [] fullLeft;
        delete [] fullRight;
      }
    }

  public:

    TestTemplate (std::istream &paramfile)
    {
      using std::cout;
      using std::endl;
      using std::numeric_limits;

      int mm, nn;

      ReadParameters (paramfile);

      D = new real_t[n];
      E = new real_t[n];
      Bidiag = new real_t[n*n];
      ReadBidiagonalMatrix ();

      Singular = new real_t[n*n];
      LASET<real_t>('a', n, n, zero, zero, Singular, n);

      if (ncLeft != n)
      {
        cout <<"Error: Invalid dimensions for left matrix." <<endl;
        exit (-1);
      }
      if (inLeft != "none")
      {
        std::ifstream infile;
				infile.open(inLeft);
        Left = Load<real_t>(mm, nn, infile);
        if (mm!=nrLeft || nn!=ncLeft)
        {
          cout <<"Error: Matrix from file '"<<inLeft
               <<"' has invalid dimensions." <<endl;
          exit (-1);
        }
        infile.close();
      }
      else if (nrLeft != 0)
      {
        // set right matrix to identity
        Left = new real_t[n*n];
        LASET<real_t>('a', n, n, zero, one, Left, n);
      }
      else
        Left = NULL;

      if (nrRight != n)
      {
        cout <<"Error: Invalid dimensions for right matrix." <<endl;
        exit (-1);
      }
      if (inRight != "none")
      {
        std::ifstream infile;
				infile.open(inRight);
        Right = Load<real_t>(mm, nn, infile);
        if (mm!=nrRight || nn!=ncRight)
        {
          cout <<"Error: Matrix from file '"<<inRight
               <<"' has invalid dimensions." <<endl;
          exit (-1);
        }
        infile.close();
      }
      else if (ncRight != 0)
      {
        // set right matrix to identity
        Right = new real_t[n*n];
        LASET<real_t>('a', n, n, zero, one, Right, n);
      }
      else
        Right = NULL;

      if (inAux != "none")
      {
        std::ifstream infile;
				infile.open(inAux);
        Auxiliary = Load<real_t>(mm, nn, infile);
        if (mm!=nrAux || nn!=ncAux)
        {
          cout <<"Error: Invalid dimensions for auxiliary matrix." <<endl;
          exit (-1);
        }
        infile.close();
      }
      else
        Auxiliary = NULL;

      info = 0;
      eps = numeric_limits<real_t>::epsilon();
    }

    virtual ~TestTemplate ()
    {
      delete [] D;
      delete [] E;
      delete [] Bidiag;
      delete [] Singular;
      delete [] Left;
      delete [] Right;
      delete [] Auxiliary;
    }

    void
    TimeIt ()
    {
      using std::cout;
      using std::endl;

      Timer timing;
      int_t ldLeft;

      if (nrLeft == 0)
        // BDSQR won't compute a left matrix, the leading dimension should be 1
        ldLeft = 1;
      else
        ldLeft = nrLeft;

      timing.Start();
      int_t ret =
        BDSQR<real_t>(uplo, n, ncRight, nrLeft, ncAux, D, E, Right, nrRight,
                      Left, ldLeft, Auxiliary, nrAux);
      timing.Stop();

      if (ret!=0)
      {
        cout <<"BDSQR failed. INFO=" <<ret <<endl;
        exit (1);
      }

      time = timing.Time();
    }

    void
    CheckResult ()
    {
      using std::cout;
      using std::endl;
      using std::ifstream;

      real_t tol(threshold * eps);

      // Check left orthogonal matrix
      if (outLeft != "none")
      {
        ifstream file;
        int_t r=nrLeft;
        int_t c=ncLeft;

        file.open(outLeft);
        if (!MatchReference (tol, file, r, c, Left, nrLeft))
        {
          cout <<"Error: computed component Q(" <<r+1 <<"," <<c+1 <<")"
               <<" and corresponding reference in file'" <<outLeft <<"'"
               <<" are not close." <<endl;
          info |= 16;
        }
        file.close();
      }
      if (nrLeft != 0)
      {
        // left matrix has been previously computed
        if (!IsOrthogonal (tol, nrLeft, ncLeft, Left, nrLeft,
                           orthogonality[0], orthogonality[1]))
        {
          cout <<"Warning: orthogonality check failed for the computed left"
               <<" matrix." <<endl;
          info |= 2;
        }
      }

      // check singular values
      if (outSingular != "none")
      {
        ifstream file;
        int_t r=n;
        int_t c=1;

        file.open(outSingular);
        if (!MatchReference (tol, file, r, c, D, 1))
        {
          cout <<"Error: computed component D(" <<r+1 <<")"
               <<" and corresponding reference in file'" <<outSingular <<"'"
               <<" are not close." <<endl;
          info |= 8;
        }
        file.close();
      }

      // check right orthogonal matrix
      if (outRight != "none")
      {
        ifstream file;
        int_t r=nrRight;
        int_t c=ncRight;

        file.open(outRight);
        if (!MatchReference (tol, file, r, c, Right, nrRight))
        {
          cout <<"Error: computed component P'(" <<r+1 <<"," <<c+1 <<")"
               <<" and corresponding reference in file'" <<outRight <<"'"
               <<" are not close." <<endl;
          info |= 32;
        }
        file.close();
      }
      if (ncRight != 0)
      {
        // right matrix has been previously computed
        if (!IsOrthogonal (tol, nrRight, ncRight, Right, nrRight,
                           orthogonality[2], orthogonality[3]))
        {
          cout <<"Warning: orthogonality check failed for the computed right"
               <<" matrix." <<endl;
          info |= 4;
        }
      }

      // check factorization
      CheckDiagonalization ();

      // check auxiliary matrix
      if (outAux != "none")
      {
        ifstream file;
        int_t r=nrAux;
        int_t c=ncAux;

        file.open(outAux);
        if (!MatchReference (tol, file, r, c, Auxiliary, nrAux))
        {
          cout <<"Error: computed component C(" <<r+1 <<"," <<c+1 <<")"
               <<" and corresponding reference in file'" <<outAux <<"'"
               <<" are not close." <<endl;
          info |=64;
        }
        file.close();
      }
    }

    void
    SaveStats ()
    {
      using std::cout;
      using std::endl;
      using std::ofstream;

      ofstream file;

      file.open(outLog);
      file <<"       start date: " <<datestart;
			file <<"       stop  date: " <<datestop;
      file <<"Elapsed time (us): " <<time <<endl;
      file <<"   ||B - Q*S*P'||= " <<residual <<endl;
      file <<"     ||Q'*Q - I||= " <<orthogonality[0] <<endl;
      file <<"     ||Q*Q' - I||= " <<orthogonality[1] <<endl;
      file <<"     ||P'*P - I||= " <<orthogonality[3] <<endl;
      file <<"     ||P*P' - I||= " <<orthogonality[2] <<endl;
      file.close();
    }

    void
    PrintStats ()
    {
      using std::cout;
      using std::endl;

			cout <<"(" <<outLog <<")" <<endl;
      cout <<"Elapsed time (us):" <<time <<endl;
      cout <<"||B - Q*S*P'||=" <<residual <<endl;
      cout <<"||Q'*Q - I||=" <<orthogonality[0] <<endl;
      cout <<"||Q*Q' - I||=" <<orthogonality[1] <<endl;
      cout <<"||P'*P - I||=" <<orthogonality[3] <<endl;
      cout <<"||P*P' - I||=" <<orthogonality[2] <<endl;
    }

		
		int_t
		Run ()
		{
			std::time_t start, stop;
			struct tm * timeinfo;

			start = std::time (NULL);
			TimeIt ();
			CheckResult ();
			stop = std::time (NULL);

			timeinfo = std::localtime (&start);
			datestart = std::asctime(timeinfo);
			timeinfo = std::localtime (&stop);
			datestop = std::asctime(timeinfo);
			if (outLog != "none")
				SaveStats ();
			
			PrintStats ();

			return info;
		}
  };

  // we need a special function to create a parametrized object
  // from the template classes of tests
  Test*
  CreateTest (std::istream &file)
  {
    using std::cout;
    using std::endl;

    std::string wprec;
    Test* test = NULL;

    file >>wprec;

    if (wprec=="double")
      test = new TestTemplate< double, double >(file);
#if 0
    if (wprec=="float")
      test = new TestTemplate< float, float >(file);
    else if (wprec=="double")
      test = new TestTemplate< double, double>(file);
    else if (wprec=="ldouble")
      test = new TestTemplate< long double, long double>(file);
    else if (wprec=="complex")
      test = new TestTemplate< float, complex<float> >(file);
    else if (wprec=="dcomplex")
      test = new TestTemplate< double, complex<double> >(file);
#ifdef REAL
    else if (wprec.compare(0, 4, "real") == 0)
      switch ( stoi(wprec.substr(4)) )
      {
      case 53:
        test = new TestTemplate< double, double>(file);
        break;
      case 128:
        test = new TestTemplate< double, double>(file);
        break;
      case 1024:
        test = new TestTemplate< double, double>(file);
        break;
      case 2048:
        test = new TestTemplate< double, double>(file);
        break;
      default :
        cout <<"Error: type '" <<wprec <<"' not supported." << endl;
        exit (-1);
      }
#endif
#ifdef MPREAL
    else if (wprec=="mpreal")
      test = new TestTemplate< mpfr::mpreal, mpfr::mpreal >(file);
#endif
    else
    {
      cout <<"Error: type '" <<wprec <<"' not supported." << endl;
      exit (-1);
    }
#endif

    return test;
  }
} // namespace TestBDSQR

//
//
//

int
main(int argc,char **argv)
{
  using std::cin;
  using std::cout;
  using std::endl;
  using std::ifstream;

  using namespace TestBDSQR;

  Test *t;

  if (argc == 1)
    t = CreateTest (cin);
  else
  {
    ifstream paramFile;
    paramFile.open(argv[1]);
    if (paramFile.fail())
    {
      cout <<"Error: Cannot open file '" <<argv[1] <<"'." <<endl;
      return -1;
    }

    t = CreateTest (paramFile);
    paramFile.close();
  }

  int errorcode = t->Run ();

  delete t;

  return errorcode;
}
