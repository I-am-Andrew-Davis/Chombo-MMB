#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IFSLICER_H_
#define _IFSLICER_H_

#include <map>

#include "Notation.H"
#include "BaseIF.H"

#include "Vector.H"
#include "REAL.H"
#include "IndexTM.H"

using std::endl;

#include "NamespaceHeader.H"

template <int dim> class IFSlicer
{
public:

  /// null constructor
  /**
      null constructor
   */
  IFSlicer();

  /// copy constructor
  /**
      copy constructor
   */
  IFSlicer(const IFSlicer<dim> & a_IFSlicer);

  /// makes a reduced IFSlicer
  /**
      makes a reduced IFSlicer
   */
  IFSlicer(IFSlicer<dim+1> const * a_slicer,
           const int             & a_fixedComp,
           const Real            & a_fixedValue);

  void setParams(IFSlicer<dim+1> const * a_slicer,
                 const int             & a_fixedComp,
                 const Real            & a_fixedValue);

  /// destructor
  /**
      destructor
   */
  virtual ~IFSlicer();

  /// Return the partial derivative evaluated at a_point.
  /**
      Represent the partial derivative as a list of integers.
      Return the partial derivative evaluated at a_point.
   */
  virtual Real value(const IndexTM<int, dim> & a_partialDerivative,
                     const IndexTM<Real,dim> & a_point)const;

  /// output function
  /**
      output function
   */
  void print(ostream& out) const;

protected:
  // The higher dimensional IFSlicer being further sliced
  IFSlicer<dim + 1> const * m_IFSlicer;

  // The component direction being held fixed and its value
  int  m_fixedComp;
  Real m_fixedValue;

private:
  //
  // These functions are made private to prevent use of the default
  // functions provided by the C++ compiler.
  //
  IFSlicer<dim>& operator=(const IFSlicer<dim>&);

};

template<> class IFSlicer<GLOBALDIM>
{
public:
  /// null constructor
  /**
      null constructor
   */
  IFSlicer();

  /// copy constructor
  /**
      copy constructor
   */
  IFSlicer(const IFSlicer<GLOBALDIM> & a_IFSlicer);

  /// if dim == GLOBALDIM, construct from BaseIF
  /**
      if dim == GLOBALDIM, construct from BaseIF
   */
  IFSlicer(const BaseIF & a_implicitFunction);

  /// destructor
  /**
      destructor
   */
  virtual ~IFSlicer();

  /// Return the partial derivative evaluated at a_point.
  /**
      Represent the partial derivative as a list of integers.
      Return the partial derivative evaluated at a_point.
   */
  virtual Real value(const IndexTM<int, GLOBALDIM> & a_partialDerivative,
                     const IndexTM<Real,GLOBALDIM> & a_point) const;

  void print(ostream& a_out) const;

protected:
  // The implicit function being sliced.  At GLOBALDIM this is just a
  // container class to start the slicing process.
  BaseIF* m_implicitFunction;

private:
  //
  // These functions are made private to prevent use of the default
  // functions provided by the C++ compiler.
  //
  IFSlicer<GLOBALDIM>& operator=(const IFSlicer<GLOBALDIM>&);
};

#include "NamespaceFooter.H"

#include "IFSlicerImplem.H"

#endif
