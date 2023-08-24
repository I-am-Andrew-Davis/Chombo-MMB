#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// SingleBlockCSAdapter.H
// Jeff J., 6 Jun 2011

#ifndef _SINGLEBLOCKCSADAPTOR_H_
#define _SINGLEBLOCKCSADAPTOR_H_

#include "MultiBlockCoordSys.H"
#include "NamespaceHeader.H"

//! \class SingleBlockCSAdapter
//! This class converts a NewCoordSys to a MultiBlockCoordSys so that
//! we can pay less attention to the differences in the two interfaces,
//! and write single- and multi-block algorithms in one framework.
class SingleBlockCSAdaptor: public MultiBlockCoordSys
{
public:

  enum Topology
  {
    TopologyPeriodic,
    TopologyDonutX
  };

  //! Creates a MultiBlockCoordSys with exactly one block, having the
  //! given coordinate system. This class assumes control of \a a_coordSys.
  SingleBlockCSAdaptor(const NewCoordSys* const a_coordSys,
                       const ProblemDomain&     a_problemDomain,
                       const RealVect&          a_domainLength,
                       const Topology&          a_topology = TopologyPeriodic);

  //! Destructor.
  virtual ~SingleBlockCSAdaptor();

  // Overridden methods.

  bool isMultiBlock() const override
    { return m_topology != TopologyPeriodic; }

  int whichBlock(const Box& a_box) const override
    { return 0; }

  int whichBlock(const IntVect& a_iv) const override
    { return 0; }

  int whichBlock(const std::string& a_name) const override
    { return 0; }

  /// Block mapping conversion function
  virtual void blockRemapping(RealVect&            a_xi_valid,
                              int&                 a_n_valid,
                              bool&                a_validExists,
                              RigidTransformation& a_extraDispl,
                              const RealVect&      a_xiSrc,
                              const IntVect&       a_iSrc,
                              int                  a_nSrc) const override;\

  ProblemDomain levelDomain() const override
    { return problemDomain(0); }

protected:

  Topology m_topology;
};


//! \class SingleBlockCSAdaptorFactory
//! Factory associated with SingleBlockCSAdaptor.
class SingleBlockCSAdaptorFactory: public MultiBlockCoordSysFactory
{
  public:

  //! Assumes control of a_coordSysFactory.
  explicit SingleBlockCSAdaptorFactory(
    const NewCoordSysFactory* const      a_coordSys,
    const RealVect                       a_domainLength = RealVect_unit,
    const SingleBlockCSAdaptor::Topology a_topology =
      SingleBlockCSAdaptor::TopologyPeriodic);

  ~SingleBlockCSAdaptorFactory();
  MultiBlockCoordSys* getCoordSys(const ProblemDomain& a_levelDomain,
                                  const RealVect& a_dx) const;

  private:

  const NewCoordSysFactory* const m_coordSysFactory;
  RealVect m_domainLength;
  SingleBlockCSAdaptor::Topology m_topology;
};

#include "NamespaceFooter.H"

#endif
