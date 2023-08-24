#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "FABView.H"
#include "parstream.H"
#include "BSplineInterp.H"
#include "BSplineVecInterp.H"
#include "CHArray.H"
#include "Box.H"
#include "RootSolver.H"
#include "CH_Timer.H"
#include <cmath>
#include <vector>
#include <algorithm>

#include "CHMatrixOps.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testBsplineInterp";
static const char* indent = "   ";
static const char* indent2 = "      ";
static bool verbose = true;

/// Compile-time constants
constexpr int c_order3 = 3;
constexpr int c_order5 = 5;

#ifdef CH_USE_DOUBLE
static Real g_precision = 1.0e-14;
#else
static Real g_precision = 1.0e-7;
#endif

///
// Parse the standard test options (-v -q -h) and
// app-specific options (-S <domain_size>) out of the command line
///
void
parseTestOptions(int argc, char* argv[])
{
  for ( int i = 1; i < argc; ++i)
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp(argv[i], "-v", 3) == 0)
            {
              verbose = true;
              //argv[i] = "" ;
            }
          else if (strncmp(argv[i] ,"-q", 3) == 0)
            {
              verbose = false;
              //argv[i] = "" ;
            }
          else if ( strncmp(argv[i], "-h", 3) == 0 )
            {
              pout() << "usage: " << pgmname << " [-hqv]" << std::endl;
              exit( 99 );
            }
        }
    }
  return;
}

int testBsplineInterp1D();
int testBsplineInterpND();

int main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions(argc, argv);
  if (verbose)
    pout () << indent2 << "Beginning " << pgmname << " ..." << std::endl;

  int status = 0;

  int testStat = testBsplineInterp1D();
  status += testStat;
  if (testStat == 0)
    {
      pout() << indent << "1D Bspline interpolation test"
             << " passed." << std::endl;
    }
  else
    {
      pout() << indent << "1D Bspline interpolation test"
           << " failed with return code " << status << endl;
    }

  testStat = testBsplineInterpND();
  status += testStat;
  if (testStat == 0)
    {
      pout() << indent << "N-dimensional Bspline interpolation test"
             << " passed." << std::endl;
    }
  else
    {
      pout() << indent << "N-dimensional Bspline interpolation test"
           << " failed with return code " << status << endl;
    }


  if (status == 0)
    {
      pout() << indent << pgmname << " passed." << std::endl;
    }
  else
    {
      pout() << indent << pgmname << " failed with return code "
             << status << std::endl;
    }

#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize();
#endif
  return status;
}


// -- interpolation testing functions --
Real cubic(Real Y)
{
  return Y*Y*Y;
}

Real Dcubic(Real Y)
{
  return 3*Y*Y;
}

Real quintic(Real Y)
{
  return Y*Y*Y*Y*Y;
}

Real Dquintic(Real Y)
{
  return 5*Y*Y*Y*Y;
}

Real sinF(Real Y)
{
  return sin(Y);
}

Real DsinF(Real Y)
{
  return cos(Y);
}

template<int Dim>
Real nLinear(stc::RVec<Dim> Y)
{ // X = product Y_i
  return Y[0];
}

template<int Dim>
Real nCubic(stc::RVec<Dim> Y)
{ // X = product Y_i^3
  stc::RVec<Dim> temp = Y*Y*Y;
  return stc::product(temp);
}

template<int Dim>
Real nCubicD(stc::RVec<Dim> Y, stc::IVec<Dim> deriv)
{ // X = product Y_i^3 // only first derivs
  Real X = 1;
  for(int i=0; i!=Dim; i++)
    {
      if(deriv[i] > 0)
        {
          Real coeff = 1.0;
          for(int c=deriv[i]; c!=0; c--)
            {
              coeff *= (Real)(3+1-c);
            }
          X *= pow(Y[i], 3-deriv[i])*coeff;
        }
      else
        {
          X *= pow(Y[i], 3);
        }
    }
  return X;
}

template<int Dim>
Real nQuintic(stc::RVec<Dim> Y)
{ // X = product Y_i^5
  stc::RVec<Dim> temp = Y*Y*Y*Y*Y;
  return stc::product(temp);
}

template<int Dim>
Real nQuinticD(stc::RVec<Dim> Y, stc::IVec<Dim> deriv)
{ // X = product Y_i^5 // only first derivs
  Real X = 1;
  for(int i=0; i!=Dim; i++)
    {
      if(deriv[i] > 0)
        {
          Real coeff = 1.0;
          for(int c=deriv[i]; c!=0; c--)
            {
              coeff *= (Real)(5+1-c);
            }
          X *= pow(Y[i], 5-deriv[i])*coeff;
        }
      else
        {
          X *= pow(Y[i], 5);
        }
    }
  return X;
}

template<int Dim>
Real nSinF(stc::RVec<Dim> Y)
{ // X = prod sin(Y_i)
  Real res = 1;
  for(int d=0; d!=Dim; d++)
    {
      res *= sin(Y[d]);
    }
  return res;
}

template<int Dim>
Real nSinDF(stc::RVec<Dim> Y, stc::IVec<Dim> deriv)
{ // X = prod sin(Y_i)
  Real res = 1;
  Real sgn;
  for(int d=0; d!=Dim; d++)
    {
      sgn = pow(-1.0, deriv[d]/2);
      if(deriv[d]%2 == 0)
        {
          res *= sgn*sin(Y[d]);
        }
      else
        {
          res *= sgn*cos(Y[d]);
        }
    }
  return res;
}

// test function for viscosity calculation from temp
Real viscFunc(Real T)
{
  return std::exp(0.746*std::log(T) +43.555/T -3259.934/(T*T) +0.13556);
}

// calculate norms
Real l2norm(Real* a_error, int size)
{
  Real val = 0;
  for(int i=0; i!=size; i++)
    {
      val += a_error[i]*a_error[i];
    }
  return sqrt(val/(Real)size);
}

Real lInfnorm(Real* a_error, int size)
{
  Real max = a_error[0]*a_error[0];
  for(int i=1; i!=size; i++)
    {
      if( a_error[i]*a_error[i] > max )
        {
          max = a_error[i]*a_error[i];
        }
    }
  return sqrt(max);
}

Real l2norm(std::vector<Real> a_error)
{
  return l2norm(a_error.data(), a_error.size());
}

Real lInfnorm(std::vector<Real> a_error)
{
  return lInfnorm(a_error.data(), a_error.size());
}

// ***********************************************************************
// test a 1D interpolation function
// ***********************************************************************
template <int Order>
int runInterpTest(const int n_pts,
                   const int n_eval,
                   const Real dx,
                   const Real xo,
                   Real(*func)(Real),
                   Real& L2,
                   Real& Linf,
                   const bool print=false,
                   const int eval_deriv=0,
                   Real(*func_deriv)(Real)=NULL)
{
  std::vector<Real> data(n_pts);
  std::vector<Real> error(n_eval);
  Real eval_dx = dx*((Real)(n_pts-1))/((Real)(n_eval-1));
  for(int i=0; i!=n_pts; i++) // set up some data
    {
      data[i] = func(i*dx + xo);
    }
  // interp
  BSplineVecInterp<Order> interp(data.data(), data.size(), xo, dx, 1);
  // solve error
  for(Real i=0; i!=n_eval; i++)
    {
      Real x = xo + i*eval_dx;
      if(eval_deriv == 0)
        {
          error[i] = func(x) - interp.interp(x);
        }
      else if(eval_deriv > 0)
        {
          error[i] = func_deriv(x) - interp.interpD(x, eval_deriv);
        }
    }
  L2 = l2norm(error);
  Linf = lInfnorm(error);
  auto maxVal = *std::max_element(data.begin(), data.end());
  Real relMax = maxVal*g_precision;
  if(print)
    {
      pout() << " - Validate order-" << Order;
      if (eval_deriv != 0)
        {
          pout() << " derivative-" << eval_deriv;
        }
      pout() << " interpolation -" << endl;
      pout() << "   (Expect floating point error relative to "
             << maxVal << ")" << endl;
      pout() << "    L_inf error: " << Linf << std::endl; 
      pout() << "    L_2   error: " << L2 << std::endl;
      pout() << std::endl;
    }
  return (Linf > relMax);
}

// ***********************************************************************
// test a 1D interpolation convergance
// ***********************************************************************
template <int Order>
int runIntpConvergTest(const int base_pts,
                       const int n_eval,
                       const Real length,
                       const Real xo,
                       Real(*func)(Real),
                       const int levels,
                       const bool print=false,
                       const int eval_deriv=0,
                       Real(*func_deriv)(Real)=NULL)
{
  std::vector<Real> L2errors(levels);
  std::vector<Real> LInferrors(levels);
  for (int L=0; L!=levels; L++)
    {
      int n_pts = (base_pts)*pow(2, L)+1;
      Real dx = length/((Real)(n_pts-1));
      runInterpTest<Order>(n_pts, n_eval, dx, xo, func,
                           L2errors[L], LInferrors[L], false,
                           eval_deriv, func_deriv);
    }
  
  if (print)
    {
      pout() << " - Order-" << Order;
      if (eval_deriv != 0)
        {
          pout() << " derivative-" << eval_deriv;
        }
      pout() << " interpolation convergence test -" << endl;
      pout() << "   N              L_2            L_inf          Order" << endl;
    }
  Real minConvg = 2*Order;
  for (int L=0; L!=levels; L++)
    {
      if (print)
        {
          pout() << "   " << (base_pts)*pow(2, L) << "   " << L2errors[L] <<
            "   " << LInferrors[L];
        }
      if (L != levels-1)
        {
          Real convRate = log(L2errors[L]/L2errors[L+1])/log(2.0);
          minConvg = std::min(minConvg, convRate);
          if (print) pout() << "   " << convRate;
        }
      if (print) pout() << endl;
    }
  // should be 1+ but using 0.5+ as a fudge factor
  return (minConvg < (0.5+Order-eval_deriv));
}

// ***********************************************************************
// test a N-Dimension interpolation function
// ***********************************************************************
template <int Dim, int Order>
int runInterpTest (const stc::IVec<Dim> n_pts,
                   const stc::IVec<Dim> n_eval,
                   const stc::RVec<Dim> dx,
                   const stc::RVec<Dim> xo,
                   Real(*func)(stc::RVec<Dim>),
                   Real& L2,
                   Real& Linf,
                   const bool print=false,
                   const stc::IVec<Dim> eval_deriv=stc::make_IVec<Dim>::zero(),
                   Real(*func_deriv)(stc::RVec<Dim>, stc::IVec<Dim>)=NULL)
{
  stc::RVec<Dim> eval_dx = dx*(n_pts-stc::make_IVec<Dim>::unit())/
    (n_eval-stc::make_IVec<Dim>::unit());
  std::vector<Real> data(stc::product(n_pts));
  std::vector<Real> error(stc::product(n_eval));
  // set data
  for(int idx=0; idx!=stc::product(n_pts); idx++)
    {
      stc::IVec<Dim> idxVec;
      int dirSz = stc::product(n_pts);
      int idxTmp = idx;
      for(int d=Dim-1; d>=0; d--)
        {
          dirSz /= n_pts[d];
          idxVec[d] = idxTmp/dirSz;
          idxTmp %= dirSz;
        }
      stc::RVec<Dim> point(xo + dx*idxVec);
      data[idx] = func(point);
    }
  
  // interpolate the data
  BSplineInterp<Dim, Order> arrInterp(
    data.data(), n_pts, 1, ColStorage, dx, xo);
  for(int idx=0; idx!=stc::product(n_eval); idx++)
    {
      stc::IVec<Dim> idxVec;
      int dirSz = stc::product(n_eval);
      int idxTmp = idx;
      for(int d=Dim-1; d>=0; d--)
        {
          dirSz /= n_eval[d];
          idxVec[d] = idxTmp/dirSz;
          idxTmp %= dirSz;
        }
      stc::RVec<Dim> point(xo + eval_dx*idxVec);
      //pout() << " at " << point << " diff " << func(point) - arrInterp.interp(point, 0) << endl;
      if(eval_deriv == stc::make_IVec<Dim>::zero()){
        error[idx] = func(point) - arrInterp.interp(point, 0);
      }
      else{
        error[idx] =
          (func_deriv(point, eval_deriv) -
           arrInterp.interpD(point, 0, eval_deriv));
      }
    }
  L2 = l2norm(error.data(), stc::product(n_eval));
  Linf = lInfnorm(error.data(), stc::product(n_eval));
  auto maxVal = *std::max_element(data.begin(), data.end());
  Real relMax = maxVal*1.0e-13;
  if(print)
    {
      pout() << " - Validate " << Dim << "-D, order-" << Order;
      if (eval_deriv != stc::make_IVec<Dim>::zero())
        {
          pout() << " derivative-" << eval_deriv;
        }
      pout() << " interpolation -" << endl;
      pout() << "   (Expect floating point error relative to "
             << maxVal << ")" << endl;
      pout() << "    L_inf error: " << Linf << std::endl;
      pout() << "    L_2   error: " << L2 << std::endl;
      pout() << std::endl;
    }
  return (Linf > relMax);
}

// ***********************************************************************
// Dump data from a N-Dimension interpolation
// ***********************************************************************
template <int Dim, int Order>
void dumpInterp (const stc::IVec<Dim> n_pts,
                 const stc::RVec<Dim> dx,
                 const stc::RVec<Dim> xo,
                 const stc::IVec<Dim> n_eval,
                 const stc::RVec<Dim> eval_dx,
                 const stc::RVec<Dim> eval_xo,
                 Real(*func)(stc::RVec<Dim>),
                 const stc::IVec<Dim> eval_deriv=stc::make_IVec<Dim>::zero())
{
  std::vector<Real> data(stc::product(n_pts));
  // set data
  for(int idx=0; idx!=stc::product(n_pts); idx++)
    {
      stc::IVec<Dim> idxVec;
      int dirSz = stc::product(n_pts);
      int idxTmp = idx;
      for(int d=Dim-1; d>=0; d--)
        {
          dirSz /= n_pts[d];
          idxVec[d] = idxTmp/dirSz;
          idxTmp %= dirSz;
        }
      stc::RVec<Dim> point(xo + dx*idxVec);
      data[idx] = func(point);
    }
  
  // interpolate the data
  BSplineInterp<Dim, Order> arrInterp(
    data.data(), n_pts, 1, ColStorage, dx, xo);
  // loop over and print
  for(int idx=0; idx!=stc::product(n_eval); idx++)
    {
      stc::IVec<Dim> idxVec;
      int dirSz = stc::product(n_eval);
      int idxTmp = idx;
      for(int d=Dim-1; d>=0; d--)
        {
          dirSz /= n_eval[d];
          idxVec[d] = idxTmp/dirSz;
          idxTmp %= dirSz;
        }
      stc::RVec<Dim> point(eval_xo + eval_dx*idxVec);
      pout() << point << "\t";
      if(eval_deriv == stc::make_IVec<Dim>::zero()){
        pout() << func(point) << "\t"
               << func(point) - arrInterp.interp(point, 0);
      }
      else{
        pout() << arrInterp.interpD(point, 0, eval_deriv);
      }
      pout() << std::endl;
    }
}

// ***********************************************************************
// test a N-D interpolation convergance
// ***********************************************************************
template <int Dim, int Order>
int runIntpConvergTest(const stc::IVec<Dim> base_pts,
                       const stc::IVec<Dim> n_eval,
                       const stc::RVec<Dim> length,
                       const stc::RVec<Dim> xo,
                       Real(*func)(stc::RVec<Dim>),
                       const int levels,
                       const bool print=false,
                       const stc::IVec<Dim> eval_deriv=stc::make_IVec<Dim>::zero(),
                       Real(*func_deriv)(stc::RVec<Dim>, stc::IVec<Dim>)=NULL)
{
  std::vector<Real> L2errors(levels);
  std::vector<Real> LInferrors(levels);
  for(int L=0; L!=levels; L++)
    {
      stc::IVec<Dim> n_pts = (base_pts)*pow(2, L) + stc::make_IVec<Dim>::unit();
      stc::RVec<Dim> n_intv = n_pts - stc::make_RVec<Dim>::unit();
      stc::RVec<Dim> dx = length/n_intv;
      runInterpTest<Dim, Order>(n_pts, n_eval, dx, xo, func,
                                L2errors[L], LInferrors[L], false,
                                eval_deriv, func_deriv);
    }

  if (print)
    {
      pout() << " - " << Dim << "-D, order-" << Order;
      if (eval_deriv != stc::make_RVec<Dim>::zero())
        {
          pout() << " derivative-" << eval_deriv;
        }
      pout() << " interpolation convergence test -" << endl;
      pout() << "   N              L_2            L_inf          Order" << endl;
    }
  Real minConvg = 2*Order;
  for (int L=0; L!=levels; L++)
    {
      if (print)
        {
          pout() << "   " << (base_pts[0])*pow(2, L) << "   " << L2errors[L] <<
            "   " << LInferrors[L];
        }
      if (L != levels-1)
        {
          Real convRate = log(L2errors[L]/L2errors[L+1])/log(2.0);
          minConvg = std::min(minConvg, convRate);
          if (print) pout() << "   " << convRate;
        }
      if (print) pout() << endl;
    }
  // should be 1+ but using 0.5+ as a fudge factor
  return (minConvg < (0.5+Order-stc::sum(eval_deriv)));
}

// ***********************************************************************
// -- Begin the tests --
// ***********************************************************************
int testBsplineInterp1D()
{
  CH_TIMERS("testBsplineInterp");
  int passed = 0;
#if CXXSTD>=14
  pout() << std::scientific << std::setprecision(6);
  Real L2, Linf;
  
  // -- 1D verification --
  if (verbose) pout() << "*** Testing bspline interpolation in 1D ***" << endl;
  {
    int n_pts = 11;
    int n_eval = 100001;
    Real dx = 0.5;
    Real xo = -1.0;
    // test re-creating cubic polynomial

    passed +=
      runInterpTest<3>(n_pts, n_eval, dx, xo, cubic, L2, Linf, verbose);
    passed +=
      runInterpTest<3>(n_pts, n_eval, dx, xo, cubic, L2, Linf, verbose, 1, Dcubic);

    // test re-creating quintic polynomial
    passed +=
      runInterpTest<5>(n_pts, n_eval, dx, xo, quintic, L2, Linf, verbose);
    passed +=
      runInterpTest<5>(n_pts, n_eval, dx, xo, quintic, L2, Linf, verbose, 1, Dquintic);
    if (verbose)
      {
        pout() << "1D verification status = " << passed << std::endl;
      }
  }
  // -- 1D convergance --
  {
    // convergance of sin
    if (verbose) pout() << "-- Convergance to sin function --" << endl;
    int levels = 5;
    int n_eval = 10001;
    Real xo = 0.0;
    Real length = 4*(atan(1) * 4);
    int base_pts = 20;
    passed +=
      runIntpConvergTest<3>(base_pts, n_eval, length, xo, sinF, levels, verbose);

    passed +=
      runIntpConvergTest<5>(base_pts, n_eval, length, xo, sinF, levels, verbose);

    passed +=
      runIntpConvergTest<3>(base_pts, n_eval, length, xo, sinF, levels, verbose, 1, DsinF);

    passed +=
      runIntpConvergTest<5>(base_pts, n_eval, length, xo, sinF, levels, verbose, 1, DsinF);

    if (verbose)
      {
        pout() << "1D convergence to sin status = " << passed << std::endl;
      }
  }
  {
    // convergance of viscosity function interpolation
    //pout() << "-- Convergance to viscous function --" << endl;
    int levels = 5;
    int n_eval = 10001;
    Real xo = 300.0;
    Real length = 700.0;
    int base_pts = 20;

    passed +=
      runIntpConvergTest<3>(base_pts, n_eval, length, xo, viscFunc, levels, verbose);

    passed +=
      runIntpConvergTest<5>(base_pts, n_eval, length, xo, viscFunc, levels, verbose);
    
    if (verbose)
      {
        pout() << "1D convergence to viscous status = " << passed << std::endl;
      }
  }
#endif
  return passed;
}

int testBsplineInterpND()
{
  CH_TIMERS("testBsplineInterp");
  int passed = 0;
#if CXXSTD>=14
  pout() << std::scientific << std::setprecision(6);

  // -- 2D verification --
  if (SpaceDim == 2)
    {
      if (verbose)
        {
          pout() << std::endl << "*** Testing bspline interpolation in 2D ***" << std::endl;
        }
    // test re-creating 2D cubic polynomial
    constexpr int dim = 2;
    stc::IVec<dim> n_pts = 21*stc::make_RVec<dim>::unit();
    stc::IVec<dim> n_eval = 101*stc::make_RVec<dim>::unit();
    stc::RVec<dim> dx = 0.2*stc::make_RVec<dim>::unit();
    stc::RVec<dim> x0 = -1.0*stc::make_RVec<dim>::unit();
    Real L2, Linf;

    passed +=
      runInterpTest<dim, c_order3> (n_pts, n_eval, dx, x0, nCubic<dim>, L2, Linf, verbose);

    passed +=
      runInterpTest<dim, c_order3> (n_pts, n_eval, dx, x0, nCubic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({1, 0}), nCubicD<dim>);

    passed +=
      runInterpTest<dim, c_order3> (n_pts, n_eval, dx, x0, nCubic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({0, 1}), nCubicD<dim>);

    passed +=
      runInterpTest<dim, c_order3> (n_pts, n_eval, dx, x0, nCubic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({1, 1}), nCubicD<dim>);

    // test re-creating 2D quintic polynomial
    passed +=
      runInterpTest<dim, c_order5> (n_pts, n_eval, dx, x0, nQuintic<dim>, L2, Linf, verbose);

    passed +=
      runInterpTest<dim, c_order5> (n_pts, n_eval, dx, x0, nQuintic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({1, 0}), nQuinticD<dim>);

    passed +=
      runInterpTest<dim, c_order5> (n_pts, n_eval, dx, x0, nQuintic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({0, 1}), nQuinticD<dim>);

    passed +=
      runInterpTest<dim, c_order5> (n_pts, n_eval, dx, x0, nQuintic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({1, 1}), nQuinticD<dim>);
    
    if (verbose)
      {
        pout() << "2D verification status = " << passed << std::endl;
      }
  }
  if (false)
    {
      // print out a bunch of values of an interpolation
      // and examine extrapolation points
      // test re-creating 2D cubic polynomial
      constexpr int dim = 2;
      stc::IVec<dim> n_pts = 11*stc::make_RVec<dim>::unit();
      stc::RVec<dim> dx = 1.0*stc::make_RVec<dim>::unit();
      stc::RVec<dim> x0 = 0.0*stc::make_RVec<dim>::unit();
      stc::IVec<dim> n_eval = (4*(10+2*6)+1)*stc::make_RVec<dim>::unit();
      stc::RVec<dim> eval_dx = 0.25*stc::make_RVec<dim>::unit();
      stc::RVec<dim> eval_x0 = -6.0*stc::make_RVec<dim>::unit();
      // stc::IVec<dim> deriv{0, 1};
    
      pout() << " - Validate Cubic interpolation -" << endl;
      pout() << "   (Expect very small floating point errors)" << endl;
      
      dumpInterp<dim, c_order5>(n_pts, dx, x0, n_eval, eval_dx, eval_x0, nLinear<dim>);
    }

  // -- 2D convergance --
  if (SpaceDim == 2)
  {
    // convergance of sin product
    if (verbose) pout() << "-- Convergance to sin function --" << endl;
    constexpr int dim = 2;
    int levels = 5;
    stc::IVec<dim> n_eval = 101*stc::make_IVec<dim>::unit();
    stc::RVec<dim> xo = stc::make_RVec<dim>::zero();
    stc::RVec<dim> length = 4*(atan(1) * 4)*stc::make_RVec<dim>::unit();
    stc::IVec<dim> base_pts = 20*stc::make_IVec<dim>::unit();
    //pout() << " - Cubic interpolation convergance test -" << endl;
    passed +=
      runIntpConvergTest<dim, c_order3> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose);
    //pout() << "   df/dx_1" << endl;
    passed +=
      runIntpConvergTest<dim, c_order3> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({1, 0}), nSinDF<dim>);
    //pout() << "   df/dx_2" << endl;
    passed +=
      runIntpConvergTest<dim, c_order3> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({0, 1}), nSinDF<dim>);
    //pout() << "   d^2f/(dx_1 d_dx_2)" << endl;
    passed +=
      runIntpConvergTest<dim, c_order3> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({1, 1}), nSinDF<dim>);

    //pout() << " - Quintic interpolation convergance test -" << endl;
    passed +=
      runIntpConvergTest<dim, c_order5> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose);
    //pout() << "   df/dx_1" << endl;
    passed +=
      runIntpConvergTest<dim, c_order5> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({1, 0}), nSinDF<dim>);
    //pout() << "   df/dx_2" << endl;
    passed +=
      runIntpConvergTest<dim, c_order5> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({0, 1}), nSinDF<dim>);
    //pout() << "   d^2f/(dx_1 d_dx_2)" << endl;
    passed +=
      runIntpConvergTest<dim, c_order5> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({1, 1}), nSinDF<dim>);
    //pout() << endl;
    if (verbose)
      {
        pout() << "2D convergence status = " << passed << std::endl;
      }
  }

  // -- 2D Inverse --
  if (SpaceDim == 2)
  {
    if (verbose) pout() << "-- Test an inverse --" << endl;
    constexpr int Dim = 2;
    stc::RVec<Dim> xo = stc::make_RVec<Dim>::zero();
    stc::RVec<Dim> length = 2*(atan(1) * 4)*stc::make_RVec<Dim>::unit();
    stc::IVec<Dim> n_pts = 20*stc::make_IVec<Dim>::unit();
    stc::RVec<Dim> n_intv = n_pts - stc::make_RVec<Dim>::unit();
    stc::RVec<Dim> dx = length/n_intv;
    int nComp = 2;
    
    std::vector<Real> data(stc::product(n_pts)*nComp);
    // set data
    for(int idx=0; idx!=stc::product(n_pts); idx++)
      {
        stc::IVec<Dim> idxVec;
        int dirSz = stc::product(n_pts);
        int idxTmp = idx;
        for(int d=Dim-1; d>=0; d--)
          {
            dirSz /= n_pts[d];
            idxVec[d] = idxTmp/dirSz;
            idxTmp %= dirSz;
          }
        stc::RVec<Dim> point(xo + dx*idxVec);
        data[idx] = point[0]; //*point[0];
        data[idx+stc::product(n_pts)] = point[1];
      }
  
    // interpolate the data using a linear function
    BSplineInterp<Dim, c_order3> arrInterp(
      data.data(), n_pts, nComp, ColStorage, dx, xo);

    // grab a test point
    stc::RVec<Dim> Xi;
    stc::RVec<Dim> X;
    Xi = 1.5*(atan(1) * 4);
    X = {arrInterp.interp(Xi, 0), arrInterp.interp(Xi, 1)};

    // Solve x given xi, when xi = f(x)
    int numIter;
    auto func =
      [&]
      (stc::RVec<Dim> xi, int comp)
        {
          return arrInterp.interp(xi, comp) - X[comp];
        };
    auto funcD =
      [&]
      (stc::RVec<Dim> xi, int comp, stc::RVec<Dim> d)
        {
          return arrInterp.interpD(xi, comp, d);
        };
    stc::RVec<Dim> Xi_inv{4.5, 5.0}; // initial guess
    Xi_inv = RootSolver::Newton(numIter, func, funcD, Xi_inv);
    //pout() << "at " << Xi_inv << " " << X << endl;
    if (stc::mag(Xi_inv - Xi) >= g_precision) passed += 1;

    stc::RVec<Dim> exact = {func(Xi_inv, 0), func(Xi_inv, 1)};
    if (stc::mag(exact) >= g_precision) passed += 1;

    if (verbose)
      {
        pout() << "2D inverse status = " << passed << std::endl;
      }
  }

  // -- 3D tests --
  if (SpaceDim >= 3)
    {
      if (verbose)
        {
          pout() << endl << endl
                 << "*** Testing bspline interpolation in SpaceDim ***" << endl;
        }
    // test re-creating cubic polynomial
    constexpr int dim = SpaceDim;
    stc::IVec<dim> n_pts = 10*stc::make_RVec<dim>::unit();
    stc::IVec<dim> n_eval = 30*stc::make_RVec<dim>::unit();
    stc::RVec<dim> dx = 0.1*stc::make_RVec<dim>::unit();
    stc::RVec<dim> x0 = 0.0*stc::make_RVec<dim>::unit();
    Real L2, Linf;

    passed +=
      runInterpTest<dim, c_order3> (n_pts, n_eval, dx, x0, nCubic<dim>, L2, Linf, verbose);

    passed +=
      runInterpTest<dim, c_order3> (n_pts, n_eval, dx, x0, nCubic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({1, 0, 0}), nCubicD<dim>);
    passed +=
      runInterpTest<dim, c_order3> (n_pts, n_eval, dx, x0, nCubic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({0, 1, 0}), nCubicD<dim>);

    passed +=
      runInterpTest<dim, c_order3> (n_pts, n_eval, dx, x0, nCubic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({0, 0, 1}), nCubicD<dim>);
    // passed +=
    //   runInterpTest<dim> (n_pts, n_eval, order, dx, x0, nCubic<dim>, L2, Linf, verbose,
    //                     stc::IVec<dim> ({1, 1, 1}), nCubicD<dim>);

    // test re-creating quintic polynomial
    passed +=
      runInterpTest<dim, c_order5> (n_pts, n_eval, dx, x0, nQuintic<dim>, L2, Linf, verbose);
    
    passed +=
      runInterpTest<dim, c_order5> (n_pts, n_eval, dx, x0, nQuintic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({1, 0, 0}), nQuinticD<dim>);

    passed +=
      runInterpTest<dim, c_order5> (n_pts, n_eval, dx, x0, nQuintic<dim>, L2, Linf, verbose,
                        stc::IVec<dim> ({0, 1, 0}), nQuinticD<dim>);

    passed +=
      runInterpTest<dim, c_order5> (n_pts, n_eval, dx, x0, nQuintic<dim>, L2, Linf, verbose,
                          stc::IVec<dim> ({0, 0, 1}), nQuinticD<dim>);
    if (verbose)
      {
        pout() << "3D verification status = " << passed << std::endl;
      }
  }
  // -- 3D convergance --
  if (SpaceDim >= 3)
  {
    // convergance of sin
    if (verbose) pout() << "-- Convergance to sin function --" << endl;
    constexpr int dim = SpaceDim;
    int levels = 2;
    stc::IVec<dim> n_eval = 31*stc::make_IVec<dim>::unit();
    stc::RVec<dim> xo = stc::make_RVec<dim>::zero();
    stc::RVec<dim> length = 4*(atan(1) * 4)*stc::make_RVec<dim>::unit();
    stc::IVec<dim> base_pts = 12*stc::make_IVec<dim>::unit();

    passed +=
      runIntpConvergTest<dim, c_order3> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose);
    passed +=
      runIntpConvergTest<dim, c_order3> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({1, 0, 0}), nSinDF<dim>);
    passed +=
      runIntpConvergTest<dim, c_order3> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({0, 1, 0}), nSinDF<dim>);
    passed +=
      runIntpConvergTest<dim, c_order3> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({0, 0, 1}), nSinDF<dim>);
    passed +=
      runIntpConvergTest<dim, c_order3> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({1, 1, 1}), nSinDF<dim>);

    passed +=
      runIntpConvergTest<dim, c_order5> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose);
    passed +=
      runIntpConvergTest<dim, c_order5> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({1, 0, 0}), nSinDF<dim>);
    passed +=
      runIntpConvergTest<dim, c_order5> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({0, 1, 0}), nSinDF<dim>);
//**FIXME want 0, 0, 1???
    passed +=
      runIntpConvergTest<dim, c_order5> (base_pts, n_eval, length, xo, nSinF<dim>, levels, verbose,
                             stc::IVec<dim> ({1, 1, 1}), nSinDF<dim>);
    if (verbose)
      {
        pout() << "3D convergence status = " << passed << std::endl;
      }
  }

// {
//   // print interpolation, using non-uniform data
//   int n_pts = 11;
//   int d2_size = 5;
//   Real data[n_pts][d2_size];
//   std::vector<Real> idp(n_pts);
//   std::vector<Real> error(n_pts);
//   Real dx = 70;
//   Real xo = 300;
//   std::pout() << "Specified points points:\n";
//   for(int i=0; i!=n_pts; i++)
//     {
//       data[i][0] = viscFunc(i*dx + xo);
//       idp[i] = i*dx + xo;
//       std::pout() << "x: " << idp[i] << "\ty: " << data[i][0] << std::endl;
//     }
  
//   BSplineVecInterp interp(data[0], n_pts, xo, dx, 5, d2_size);
//   std::pout() << "Interpolation points:\n";
//   for(Real i=0; i<=(dx*(n_pts-1)); i+=0.25*dx)
//     {
//       Real x = i+xo;
//       std::pout() << "x: " << x << "  \ty: " << interp.interp(x)
//                 << "\tTrue: " << viscFunc(x) << "\terror: "
//                 << viscFunc(x) - interp.interp(x) << std::endl;
//     }
// }
  
// {
//   std::pout() << std::endl << std::endl;
//   //std::vector<Real> data = {0, 0.5, sqrt(3.0)/2.0, 1, sqrt(3.0)/2.0, 0.5, 0};
//   //Real dx = 3.14159265/data.size();
//   std::vector<Real> data(7);
//   Real dx = 1.0;
//   for(int i=0; i!=data.size(); i++)
//     {
//       data[i] = pow(i*dx, 3.0);
//       std::pout() << data[i] << " ";
//     }
//   std::pout() << std::endl;
//   BSplineVecInterp interp(data.data(), data.size(), 0.0, dx, 3, 1);
//   int deriv = 1;
//   // solve error
//   for(Real i=0; i!=data.size()*2-1; i++)
//     {
//       Real point = i*dx*0.5;
//       std::pout() << "interp deriv: " << interp.interpD(point, deriv)
//                 << " at " << point << std::endl;
//     }
//   std::pout() << std::endl << std::endl;
// }
  
  
// {
//   std::pout() << std::setprecision(2) << std::endl << std::endl;
//   // testing with CH_arrays
//   const int dim = 2;
//   int size_x = 7;
//   int size_y = 7;
//   int comp = 1;
//   // Real tDat[size_x*size_y] = {1,1,1,1,1,1,1,
//   //                             1,1,1,1,1,1,1,
//   //                             1,1,1,1,1,1,1,
//   //                             1,1,1,2,1,1,1,
//   //                             1,1,1,1,1,1,1,
//   //                             1,1,1,1,1,1,1,
//   //                             1,1,1,1,1,1,1};
//   // Real tDat[size_x*size_y] = {0, 0.5, sqrt(3.0)/2.0, 1, sqrt(3.0)/2.0, 0.5, 0,
//   //                             0, 0.5, sqrt(3.0)/2.0, 1, sqrt(3.0)/2.0, 0.5, 0,
//   //                             0, 0.5, sqrt(3.0)/2.0, 1, sqrt(3.0)/2.0, 0.5, 0,
//   //                             0, 0.5, sqrt(3.0)/2.0, 1, sqrt(3.0)/2.0, 0.5, 0,
//   //                             0, 0.5, sqrt(3.0)/2.0, 1, sqrt(3.0)/2.0, 0.5, 0,
//   //                             0, 0.5, sqrt(3.0)/2.0, 1, sqrt(3.0)/2.0, 0.5, 0,
//   //                             0, 0.5, sqrt(3.0)/2.0, 1, sqrt(3.0)/2.0, 0.5, 0,};
//   Real tDat[size_x*size_y] = {0, 0, 0, 0, 0, 0, 0,
//                               1, 1, 1, 1, 1, 1, 1,
//                               8, 8, 8, 8, 8, 8, 8,
//                               27, 27, 27, 27, 27, 27, 27,
//                               64, 64, 64, 64, 64, 64, 64,
//                               125, 125, 125, 125, 125, 125, 125,
//                               216, 216, 216, 216, 216, 216, 216};
//   // Real tDat[size_x*size_y] = {0, 0, 0, 0, 0, 0, 0,
//   //                             0, 1, 8, 27, 64, 125, 216,
//   //                             0, 8, 64, 216, 512, 1000, 1728,
//   //                             0, 27, 216, 729, 1728, 3375, 5832,
//   //                             0, 64, 512, 1728, 4096, 8000, 13824,
//   //                             0, 125, 1000, 3375, 8000, 15625, 27000,
//   //                             0, 216, 1728, 5832, 13824, 27000, 46656};
//   // Real tDat[size_x*size_y*comp] = {1,1,1,1,1,1,1,
//   //                                  1,1,1,1,1,1,1,
//   //                                  1,1,1,1,1,1,1,
//   //                                  1,1,1,1,1,1,1,
//   //                                  1,1,1,1,1,1,1,
//   //                                  1,1,1,1,1,1,1,
//   //                                  1,1,1,1,1,1,1,
//   //                                  1,2,3,4,5,6,7,
//   //                                  3,4,5,6,7,8,9,
//   //                                  5,6,7,8,9,10,11,
//   //                                  7,8,9,10,11,12,13,
//   //                                  9,10,11,12,13,14,15,
//   //                                  11,12,13,14,15,16,17,
//   //                                  13,14,15,16,17,18,19};
//   // Real tDat[size_x*size_y] = {1,2,3,4,5,6,7,
//   //                             3,4,5,6,7,8,9,
//   //                             5,6,7,8,9,10,11,
//   //                             7,8,9,10,11,12,13,
//   //                             9,10,11,12,13,14,15,
//   //                             11,12,13,14,15,16,17,
//   //                             13,14,15,16,17,18,19};
//   // print
//   // for(int c=0; c!=comp; c++)
//   //   {
//   //     for(int y=0; y!=size_y; y++)
//   //       {
//   //         for(int x=0; x!=size_x; x++)
//   //           {
//   //             std::pout() << tDat[x+size_x*y+c*size_x*size_y] << ", ";
//   //           }
//   //         std::pout() << std::endl;
//   //       }
//   //     std::pout() << std::endl;
//   //   }
//   // defualt array
//   IntVect size(size_x, size_y);
//   CHArray<Real, dim+1, ArZeroCol> tArr;
//   tArr.define(size, comp);
//   for(int c=0; c!=comp; c++)
//     {
//       for(int y=0; y!=size_y; y++)
//         {
//           for(int x=0; x!=size_x; x++)
//             {
//               IntVect pointI(x, y);
//               tArr(pointI, c) = tDat[x+y*size_x+c*size_x*size_y];
//               std::pout() << tArr(pointI, c) << ", ";
//             }
//           std::pout() << std::endl;
//         }
//       std::pout() << std::endl;
//     }
//   // std::pout() << tArr << std::endl;
//   stc::IVec<dim> sizeV({size_x, size_y});
//   stc::RVec<dim> dx({1.0, 1.0});//dx({3.14159265/size_x, 3.14159265/size_y});
//   BSplineInterp<dim> arrInterp(tArr.begin(), sizeV, comp, ColStorage, 3, dx);
//   //BSplineInterp<dim> arrInterp(tDat, sizeV, comp);
    
//   // stc::RVec<dim> point({4.9,2.3});
//   stc::RVec<dim> deriv({0,1});
//   // std::pout() << "interp: " << arrInterp.interpD(point,1,deriv) << " at " << point << std::endl;
//   for(int c=0; c!=comp; c++)
//     {
//       for(int y=0; y!=size_y*2-1; y++)
//         {
//           for(int x=0; x!=size_x*2-1; x++)
//             {
//               stc::RVec<dim> point({0.5*dx[0]*x, 0.5*dx[1]*y});
//               std::pout() << "interp deriv: " << arrInterp.interpD(point, c, deriv)
//                         << " at " << point << std::endl;
//             }
//           std::pout() << std::endl;
//         }
//       std::pout() << std::endl << std::endl;
//     }
// }
#endif
  return passed;
}
