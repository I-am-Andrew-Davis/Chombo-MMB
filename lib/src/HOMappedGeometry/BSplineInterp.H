#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___ 
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/******************************************************************************/
/**
 * \file BSplineInterp.cpp
 *
 * \brief Class for interpolation between points in arrays
 *
 *//*+*************************************************************************/

#ifndef _BSPLINEINTERP_H_
#define _BSPLINEINTERP_H_

#include "BSplineVecInterp.H"

#define USE_STCVECTOR
#ifdef USE_STCVECTOR
#include "StcVector.H"
#endif

#include "NamespaceHeader.H"

// array storage options
enum storageType
{
  ColStorage = 0, // first dir is continuous
  RowStorage = 1, // last is continuous
};

// data needed for a multi dimension array
template<int Dim=SpaceDim>
struct MDarray
{
  Real* data = nullptr;
  stc::IVec<Dim> size;
  stc::IVec<Dim> stride;
  int nComp;

  int elements() const
    {
      return nComp*product(size);
    }
  
  int linearSize() const
    {
      return elements()*sizeof(Real);
    }
  
  void linearOut(void* const a_outBuf) const
    {
      CH_assert(Misc::isAligned<Real>(a_outBuf));
      Real* buf = (Real*)a_outBuf;
      for (int i=0; i!=nComp*product(size); ++i)
        {
          *buf = data[i];
          ++buf;
        }
      //a_outBuf = (void* const)this->data;
    }
  
  void linearIn(const void* a_inBuf) const
    {
      CH_assert(Misc::isAligned<Real>(a_inBuf));
      Real* buf = (Real*)a_inBuf;
      for (int i=0; i!=nComp*product(size); ++i)
        {
          data[i] = *buf;
          ++buf;
        }
      //this->data = (Real*)a_inBuf;
    }
};

// A type for breaking MD array into 1D slices,
// for reduction into set of 1D interps
struct ArrSlice
{
  Real* data = nullptr;
  int dir;
  int stride;
  int size;
};


/*******************************************************************************
 */
///  Multi-dimensional b-spline interpolation of array data
///  Components are interpolated seperatly
///  Data must be uniformly spaced in current implementation
///  Supports use of linear(1), cubic(3), and quintic(5) interpolation
///  Solves for continious fit through first Order-1 derivitives
///  providing for O( x^(Order+1) )
///  Requires minimum of 6 points for cubic, 7 for quintic
/**
 *   // Example usage : see testBsplineInterp
 *
 *  \tparam[Order]      Polynomial order (must be odd, default 5)
 *
 */
/**
 ******************************************************************************/
template<int Dim, int Order>
class BSplineInterp
{
  // name vectors in spline space
  using ISVec = stc::IVec<Dim>;  // int spline vect
  using RSVec = stc::RVec<Dim>;  // real spline vect
  
  
/*==============================================================================
 * Data Members
 *============================================================================*/

protected:

  // Constants
  static constexpr ISVec c_esize{
    1, Order+1, Order+1, Order+1, Order+1, Order+1 };
  static constexpr ISVec c_estride = stc::cstride(c_esize);
  static constexpr int c_numS = stc::product(c_esize);

  MDarray<Dim> m_interpArr;              ///< pointer to interpolation data
  MDarray<Dim> m_ctrlArr;                ///< store control points
  MDarray<Dim> cpEval;                   ///< for evaluation
  int m_sIdx[Dim][stc::product(c_esize)];
  
  int m_arrStorage;
  std::vector<std::vector<Real> > m_T;   ///< set of knot vectors
  std::vector<Real> m_dT;                ///< uniform knot vectors spacing
  RSVec m_dx;                            ///< physical spacing
  RSVec m_x0;                            ///< physical start location
  int m_nComp;                           ///< components for interp
  int m_numProc;
  int m_startProc;
  
/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:
  /// defult constructor
  BSplineInterp()
    {
    };

  /*--------------------------------------------------------------------*/
  //  Constructor
  /** \param[in]  a_interpArr   pointer to interpolate data array
   *  \param[in]  a_size        vector of array size
   *  \param[in]  a_numComp     number of components, each done seperatly
   *  \param[in]  a_arrStorage  row or column
   *  \param[in]  a_dx          physical spacing between interp points
   *  \param[in]  a_xstart      physical coordinates of first interp point
   *  \param[in]  a_numProc     Number of processors for parallel construction
   *  \param[in]  a_startProc   First processor for parallel construction
   *//*-----------------------------------------------------------------*/
  BSplineInterp(Real*                         a_interpArr,
                const stc::Vector<int, Dim>&  a_size,
                const int                     a_numComp,
                const int                     a_arrStorage = ColStorage,
                const stc::Vector<Real, Dim>& a_dx =
                  stc::make_RVec<Dim>::unit(),
                const stc::Vector<Real, Dim>& a_xstart =
                  stc::make_RVec<Dim>::zero(),
                const int                     a_numProc = 0,
                const int                     a_startProc = 0)
    {
      define(a_interpArr,
             a_size,
             a_numComp,
             a_arrStorage,
             a_dx,
             a_xstart,
             a_numProc,
             a_startProc);
    };
  
  /*--------------------------------------------------------------------*/
  //  Weak constructor
  /** \param[in]  a_interpArr   pointer to interpolate data array
   *  \param[in]  a_size        vector of array size
   *  \param[in]  a_numComp     number of components, each done seperatly
   *  \param[in]  a_arrStorage  row or column
   *  \param[in]  a_dx          physical spacing between interp points
   *  \param[in]  a_xstart      physical coordinates of first interp point
   *//*-----------------------------------------------------------------*/
  void define(Real*                         a_interpArr,
              const stc::Vector<int, Dim>&  a_size,
              const int                     a_numComp,
              const int                     a_arrStorage = ColStorage,
              const stc::Vector<Real, Dim>& a_dx =
                stc::make_RVec<Dim>::unit(),
              const stc::Vector<Real, Dim>& a_xstart =
                stc::make_RVec<Dim>::zero(),
              const int                     a_numProc = 0,
              const int                     a_startProc = 0)
    {
      m_arrStorage = a_arrStorage;
      m_dx = a_dx;
      m_x0 = a_xstart;
      m_nComp = a_numComp;
      // set the interp array
      m_interpArr.data = a_interpArr;
      m_interpArr.size = a_size;
      m_interpArr.nComp = a_numComp;
      // parallel construction options
      m_numProc = a_numProc;
      m_startProc = a_startProc;
      // solve for the strides
      sizeToStride(m_interpArr.size, m_interpArr.stride, a_arrStorage);
      // size the control point array
      allocArr(m_ctrlArr,
               m_interpArr.size + (Order-1)*stc::make_IVec<Dim>::unit(),
               m_nComp);
      // allocate space for interpolation
      allocArr(cpEval, c_esize, 1);
      // create knot vectors
      knots();
      // solve the control points
      solveAllCP();
      // Define spline indices
      for (int s = 0; s != c_numS; ++s)
        {
          int sIdx = 0;
          int idx = s;
          for (int d = 0; d != Dim; ++d)
            {
              sIdx += (idx%c_esize[d])*m_ctrlArr.stride[d];
              idx /= c_esize[d];
            }
          m_sIdx[0][s] = sIdx;
        }
      for (int d = 1; d != Dim; ++d)
        {
          for (int s = 0; s != c_numS; ++s)
            {
              int sIdx = 0;
              int idx = s;
              for (int d2 = 0; d2 != Dim; ++d2)
                {
                  if (d2 != d)
                    {
                      sIdx += (idx%c_esize[d2])*cpEval.stride[d2];
                      idx /= c_esize[d2];
                    }
                }
              m_sIdx[d][s] = sIdx;
            }
        }
    };

  /// Copy constructor disabled
  BSplineInterp(const BSplineInterp&) = delete;

  /// Move constructor disabled
  BSplineInterp(BSplineInterp&&) = delete;

  /// Assignment constructor disabled
  BSplineInterp& operator=(const BSplineInterp&) = delete;

  /// Move assignment constructor disabled
  BSplineInterp& operator=(BSplineInterp&&) = delete;

  /*--------------------------------------------------------------------*/
  //  Destructor
  /** 
   *//*-----------------------------------------------------------------*/
  ~BSplineInterp()
    {
      // de-allocate the arrays
      deallocArr(m_ctrlArr);
      deallocArr(cpEval);
    };

/*==============================================================================
 * Public member functions
 *============================================================================*/
  
  /*--------------------------------------------------------------------*/
  //  Interpolate point in domain
  /** \param[in]  a_xi     location to interp
   *                       must be within the size of the baseFab box
   *  \return              interpolated value
   *//*-----------------------------------------------------------------*/
  Real interp(const stc::Vector<Real, Dim>& a_xi, const int a_comp) const
    {
      return interpSpline(a_xi, a_comp);
    }
      
  /*--------------------------------------------------------------------*/
  //  Interpolate point in domain by evaluating the spline
  /** \param[in]  a_xi     location to interp
   *                       must be within the size of the baseFab box
   *  \return              interpolated value
   *//*-----------------------------------------------------------------*/
  Real interpSpline(const stc::Vector<Real, Dim>& a_xi, const int a_comp) const
    {
      // CH_TIME("BSplineInterp::interp");
      CH_TIMERS("BSplineInterp::interp");
      CH_TIMER("BSplineInterp::interp::1", tm1);
      CH_TIMER("BSplineInterp::interp::2", tm2);
      CH_TIMER("BSplineInterp::interp::3", tm3);
      CH_TIMER("BSplineInterp::interp::4", tm4);
      // // check interp point is in bounds
      // CH_assert(a_xi >= m_x0 - 0.00000001);
      // CH_assert(a_xi <= m_x0 + (m_size-1)*m_dx + 0.00000001);

      CH_START(tm1);
      // normalize point
      const RSVec t = (a_xi - m_x0)/m_dx;
      // get knot indcies
      ISVec Tidx;
      for (int d = 0; d != Dim; ++d)
        {
          Real tt = t[d];
          // Check to make sure in bounds
          // If not in bounds, extrapolate using the nearest spline segment
          if (tt < 0.) tt = 0;
          else if (tt >= m_interpArr.size[d]-1) tt = m_interpArr.size[d]-2;
          // get the index
          Tidx[d] = getKnotIndx(tt, m_T[d], m_dT[d]);
        }
      CH_STOP(tm1);
      CH_START(tm2);

      // Solve basis functions
      stc::RVec<Order+1> bf;
      BSpline::BasisFuncs<Order, 0>::eval(
        bf.dataPtr(), Tidx[0], t[0], m_T[0]);
      
      // Solve the first direction seperate because it pulls from ctrlArr
      const int ctrlStride = m_ctrlArr.stride[0];
      int sIdx_base = stc::sum(m_ctrlArr.stride*(Tidx - Order));
      sIdx_base += a_comp*stc::product(m_ctrlArr.size);
      for (int s = 0; s != c_numS; ++s)
        {
          const int sIdx = m_sIdx[0][s] + sIdx_base;
          cpEval.data[s] = BSpline::evalPoint<Order>(&m_ctrlArr.data[sIdx],
                                                     ctrlStride,
                                                     bf.dataPtr());
        }
      CH_STOP(tm2);
      CH_START(tm3);

      // now other directions "recursively"
      int slcPerDir = c_numS;
      for (int d = 1; d != Dim; ++d)
        {
          // Solve basis functions
          BSpline::BasisFuncs<Order, 0>::eval(
            bf.dataPtr(), Tidx[d], t[d], m_T[d]);

          // solve over each 1D slice
          slcPerDir /= (Order + 1);
          const int cpStride = c_estride[d];
          for (int s=0; s!=slcPerDir; s++)
            {
              const int sIdx = m_sIdx[d][s];
              Real* data = cpEval.data + sIdx;
              data[0] = BSpline::evalPoint<Order>(data,
                                                  cpStride,
                                                  bf.dataPtr());
            }
        }
      CH_STOP(tm3);
      CH_START(tm4);
      //CH_assert(isfinite(cpEval.data[0]));
      if (!isfinite(cpEval.data[0]))
        {
          MayDay::Warning("BSpline evaluated is not finite");
        }
      CH_STOP(tm4);
      return cpEval.data[0];
    };


  /*--------------------------------------------------------------------*/
  //  derivitive of interpolated point in domain
  /** \tparam Deriv...     Derivatives to find in order of dimensions
   *  \param[in]  a_xi     location to interp
   *                       must be within the size of the baseFab box
   *  \param[in]  a_deriv  which derivitive to find
   *                       only valid for derivitives < spline order
   *  \return              interpolated value
   *//*-----------------------------------------------------------------*/
  Real interpD(const stc::Vector<Real, Dim>& a_xi,
               const int                     a_comp,
               const stc::Vector<int, Dim>   a_deriv) const
    {
      // CH_TIME("BSplineInterp::interpD");
      CH_TIMERS("BSplineInterp::interpD");
      CH_TIMER("BSplineInterp::interpD::1", tm1);
      CH_TIMER("BSplineInterp::interpD::2", tm2);
      CH_TIMER("BSplineInterp::interpD::3", tm3);
      CH_TIMER("BSplineInterp::interpD::4", tm4);
      // // check interp point is in bounds
      // CH_assert(a_xi >= m_x0 - 0.00000001);
      // CH_assert(a_xi <= m_x0 + (m_size-1)*m_dx + 0.00000001);

      CH_START(tm1);
      // normalize point
      RSVec t = (a_xi-m_x0)/m_dx;
      // get knot indcies
      ISVec Tidx;
      for (int d=0; d!=Dim; d++)
        {
          Real tt = t[d];
          // Check to make sure in bounds
          // If not in bounds, extrapolate using the nearest spline segment
          if (tt < 0) tt = 0;
          else if (tt >= m_interpArr.size[d]-1) tt = m_interpArr.size[d]-2;
          // get the index
          Tidx[d] = getKnotIndx(tt, m_T[d], m_dT[d]);
        }
      CH_STOP(tm1);
      CH_START(tm2);

      // Solve basis functions
      stc::RVec<Order+1> bf;
      switch(a_deriv[0])
        {
        case 0:
          BSpline::BasisFuncs<Order, 0>::eval(
            bf.dataPtr(), Tidx[0], t[0], m_T[0]);
          break;
        case 1:
          BSpline::BasisFuncs<Order, 1>::eval(
            bf.dataPtr(), Tidx[0], t[0], m_T[0]);
          break;
        default:
          BSpline::BasisFuncs<Order>::eval(
            bf.dataPtr(), Tidx[0], t[0], m_T[0], a_deriv[0]);
          break;
        }

      // Solve the first direction seperate because it pulls from ctrlArr
      const int ctrlStride = m_ctrlArr.stride[0];
      int sIdx_base = stc::sum(m_ctrlArr.stride*(Tidx - Order));
      sIdx_base += a_comp*stc::product(m_ctrlArr.size);
      for (int s = 0; s != c_numS; ++s)
        {
          const int sIdx = m_sIdx[0][s] + sIdx_base;
          switch(a_deriv[0])
            {
            case 0:
              cpEval.data[s] = BSpline::PointDeriv<Order, 0>::eval(
                Tidx[0],
                m_T[0],
                m_dx[0],
                &m_ctrlArr.data[sIdx],
                ctrlStride,
                bf.dataPtr());
              break;
            case 1:
              cpEval.data[s] = BSpline::PointDeriv<Order, 1>::eval(
                Tidx[0],
                m_T[0],
                m_dx[0],
                &m_ctrlArr.data[sIdx],
                ctrlStride,
                bf.dataPtr());
              break;
            default:
              cpEval.data[s] = BSpline::PointDeriv<Order>::eval(
                Tidx[0],
                m_T[0],
                m_dx[0],
                &m_ctrlArr.data[sIdx],
                ctrlStride,
                bf.dataPtr(),
                a_deriv[0]);
              break;
            }
        }
      CH_STOP(tm2);
      CH_START(tm3);

      // now other directions "recursively"
      int slcPerDir = c_numS;
      for (int d = 1; d != Dim; ++d)
        {
          slcPerDir /= (Order + 1);
          const int cpStride = c_estride[d];
          // Solve basis functions
          switch (a_deriv[d])
            {
            case 0:
              BSpline::BasisFuncs<Order, 0>::eval(
                bf.dataPtr(), Tidx[d], t[d], m_T[d]);
              // solve over each 1D slice
              // slc is a counter, not actual index to data
              for (int s = 0; s != slcPerDir; ++s)
                {
                  // cut the slices from each MDarray
                  const int sIdx = m_sIdx[d][s];
                  Real* data = cpEval.data + sIdx;
                  data[0] = BSpline::PointDeriv<Order, 0>::eval(
                    Tidx[d],
                    m_T[d],
                    m_dx[d],
                    data,
                    cpStride,
                    bf.dataPtr());
                }
              break;
            case 1:
              BSpline::BasisFuncs<Order, 1>::eval(
                bf.dataPtr(), Tidx[d], t[d], m_T[d]);
              // solve over each 1D slice
              // slc is a counter, not actual index to data
              for (int s = 0; s != slcPerDir; ++s)
                {
                  // cut the slices from each MDarray
                  const int sIdx = m_sIdx[d][s];
                  Real* data = cpEval.data + sIdx;
                  data[0] = BSpline::PointDeriv<Order, 1>::eval(
                    Tidx[d],
                    m_T[d],
                    m_dx[d],
                    data,
                    cpStride,
                    bf.dataPtr());
                }
              break;
            default:
              BSpline::BasisFuncs<Order>::eval(
                bf.dataPtr(), Tidx[d], t[d], m_T[d], a_deriv[d]);
              // solve over each 1D slice
              // slc is a counter, not actual index to data
              for (int s = 0; s != slcPerDir; ++s)
                {
                  // cut the slices from each MDarray
                  const int sIdx = m_sIdx[d][s];
                  Real* data = cpEval.data + sIdx;
                  data[0] = BSpline::PointDeriv<Order>::eval(
                    Tidx[d],
                    m_T[d],
                    m_dx[d],
                    data,
                    cpStride,
                    bf.dataPtr(),
                    a_deriv[d]);
                }
              break;
            }
        }
      CH_STOP(tm3);
      CH_START(tm4);
      //CH_assert(isfinite(cpEval.data[0]));
      if (!isfinite(cpEval.data[0]))
        {
          MayDay::Warning("BSpline derivative evaluated is not finite");
        }
      CH_STOP(tm4);
      return cpEval.data[0];
    }

/*==============================================================================
 * Protected member functions
 *============================================================================*/
protected:

  /*--------------------------------------------------------------------*/
  //   allocate an MDarray on the heap, set values
  /**
   *//*-----------------------------------------------------------------*/
  void allocArr(MDarray<Dim>& arr, const ISVec& size, const int nComp)
    {
      // make sure to delete prior data
      deallocArr(arr);
      // now set all the stuff
      arr.size = size;
      arr.nComp = nComp;
      sizeToStride(arr.size, arr.stride, m_arrStorage);
      arr.data = new Real[arr.elements()];
      for (int i=0; i!=arr.elements(); i++)
        {
          arr.data[i] = 0.0;
        }
    }

  /*--------------------------------------------------------------------*/
  //  deallocate an MDarray, clear mem
  /**
   *//*-----------------------------------------------------------------*/
  void deallocArr(MDarray<Dim>& arr)
    {
      if (arr.data != nullptr)
        {
          delete [] arr.data;
          // set all to zero
          // arr.data = nullptr;
          arr.size = stc::make_IVec<Dim>::zero();
          arr.stride = stc::make_IVec<Dim>::zero();
          arr.nComp = 0;
        }
    }
  
  /*--------------------------------------------------------------------*/
  //  Create the set of knot vectors, 1 per direction
  /** Each vector is uniformly spaced of unit stride
   *  There are repeated knots at the end points as needed
   *  If a direction is of size 1, then only single knot at 0
   *//*-----------------------------------------------------------------*/
  void knots()
    {
      // set number of knot vects, then set each vector along FAB directions
      m_T.resize(Dim);
      m_dT.resize(Dim);
      for (int i=0; i!=Dim; i++)
        {
          // set knot length
          m_T[i].resize(m_interpArr.size[i] + 2*Order);
          m_dT[i] = 1.0;
          // fill vector with uniform points
          for (int j=0; j!=Order; j++)
            {
              m_T[i][j] = 0.0;
            }
          for (int j=Order; j!=m_interpArr.size[i]+Order; j++)
            {
              m_T[i][j] = m_dT[i]*(j-Order);
            }
          for (int j=m_interpArr.size[i]+Order; j!=m_T[i].size(); j++)
            {
              m_T[i][j] = m_dT[i]*(m_interpArr.size[i] - 1);
            }
        }
    }

  /// get index in knot vector for any t
  int getKnotIndx(const Real t,
                  const std::vector<Real>& T,
                  const Real dT) const
    {
      return static_cast<int>((t - std::fmod(t,dT))/dT) + Order;
    }
  
  /// solve control net over all interpolation points
  /*--------------------------------------------------------------------*/
  //  Solve the MDarray of control points
  /** This is done as a recursive set of 1D problems
   *  FIXME - This needs to truly work in MPI
   *          The current quick fix is to solve on one proc and distribute
   *//*-----------------------------------------------------------------*/
  void solveAllCP()
    {
      CH_TIME("BSplineInterp::solveAllCP");
#ifdef CH_MPI
      if ((procID() == m_startProc))
#endif
        {
          // pointers to correct data to use at given point
          // these are single component
          MDarray<Dim> interp, contrl;
          // create the slices
          ArrSlice interpSl, contrlSl;
          // storage for basis functions
          stc::RVec<Order+1> bf;
      
          // solve for each component seperatly
          for (int comp=0; comp!=m_interpArr.nComp; comp++)
            {
              // set the interpolation data pointer, shift pointer by component
              interp = m_interpArr;
              interp.data = &interp.data[product(interp.size)*comp];
              // set the control point pointer, shift pointer by component
              contrl = m_ctrlArr;
              contrl.data = &contrl.data[product(contrl.size)*comp];
              // solve for CP in each dimension
              // for (int dir=Dim-1; dir>=0; dir--)
              for (int dir=0; dir!=Dim; dir++) 
                {
                  // solve over each 1D slice
                  // loop over each face normal to the slice direction
                  ISVec loopRegion(interp.size - stc::make_IVec<Dim>::unit());
                  loopRegion[dir] = 0;
                  stc::nestedLoop(
                    stc::make_IVec<Dim>::zero(), loopRegion,
                    [&]
                    (const ISVec& iv)
                    {
                      // cut the slices from each MDarray
                      interpSl.data = &interp.data[getArrIndex(interp.stride, iv, 0)];
                      interpSl.stride = interp.stride[dir];
                      interpSl.size = interp.size[dir];
                      interpSl.dir = dir;

                      contrlSl.data = &contrl.data[getArrIndex(contrl.stride, iv, 0)];
                      contrlSl.stride = contrl.stride[dir];
                      contrlSl.size = contrl.size[dir];
                      contrlSl.dir = dir;
                      
                      // solve a single slice, fill the control points
                      BSpline::solveCP<Order>(bf.dataPtr(),
                                              interpSl.data,
                                              contrlSl.data,
                                              interpSl.size,
                                              Order,
                                              m_dT[dir],
                                              m_T[dir],
                                              interpSl.stride,
                                              contrlSl.stride);
                    });
                  interp.data = contrl.data;
                  interp.stride = contrl.stride;
                  interp.size[dir] = contrl.size[dir];
                }
            }
        }
      //broadcast cp to all procs
      broadcast(m_ctrlArr, m_startProc);
    }
 
  /*--------------------------------------------------------------------*/
  ///   calculate array index
  /**
   *//*-----------------------------------------------------------------*/
  int getArrIndex(const ISVec a_stride,
                  const ISVec a_idx,
                  const int a_comp) const
    {
      CH_TIME("BSplineInterp::getArrIndex");
      // should fix depending on storageType
      int idx = 0;
      for (int i = 0; i!=Dim; i++)
        {
          idx += a_idx[i]*a_stride[i];
        }
      //idx += a_comp*product(a_stride);
      return idx;
    }

  void getMDIndex(ISVec& a_idxVec,
                  const int a_index,
                  const ISVec& a_stride)
    {
      CH_assert(m_arrStorage == ColStorage); // FIXME
      int idxTmp = a_index;
      for (int d=Dim-1; d>=0; d--)
        {
          a_idxVec[d] = idxTmp/a_stride[d];
          idxTmp %= a_stride[d];
        }
    }
  /*--------------------------------------------------------------------*/
  /// calculate stride from the size
  /**
   *//*-----------------------------------------------------------------*/
  void sizeToStride(const ISVec a_size,
                    ISVec& a_stride,
                    const int a_arrStorage = ColStorage) const
    {
      switch (a_arrStorage)
        {
        case ColStorage:
          a_stride = stc::stride(a_size);
          return;
        case RowStorage:
          a_stride = stc::stride(ISVec(a_size).reverse());
          return;
        }
    }

  void incrSliceIndex(int& a_index,
                      const ISVec& array_size,
                      const ISVec& array_stride,
                      const ISVec& loop_size,
                      const int dir)
    {
      CH_assert(loop_size <= array_size);
      ISVec loopRegion(loop_size);
      loopRegion[dir] = 1;
      // get a MD index
      ISVec indxVec;
      getMDIndex(indxVec, a_index, array_stride);
      // call dimension recursive increment
      incrMDVec(indxVec, loopRegion, 0);
      // increment to start of next slice
      a_index = getArrIndex(array_stride, indxVec, 0);
      pout() << indxVec << " " << a_index << std::endl;
    }

  bool incrMDVec(ISVec& idx, const ISVec& region, const int d)
    {
      // some checks
      //CH_assert((d < Dim)&&(d > 0));
      CH_assert(idx[d] <= region[d]);
      CH_assert(idx[d] >= 0);
      CH_assert(region[d] > 0);
      // increment
      idx[d] += 1;
      // check to see if that pushed over a bound, and correct
      if (idx[d] == region[d])
        {
          idx[d] = 0;
          incrMDVec(idx, region, d+1);
        }
      return true;
    }
  
  /*--------------------------------------------------------------------*/
  //  Take a 1D slice, spans an entire direction length
  //  Should only be used from solveAllCP
  /** \param[in]  a_base  multi-dimensional array to take slice from
   *  \param[in]  a_dir   direction of slice
   *  \param[in]  a_index data index for for slice to start at
   *  \param[in]  a_comp  the component to take the slice from
   *  \return             Resulting array slice
   *//*-----------------------------------------------------------------*/
  ArrSlice sliceArr(const MDarray<Dim> a_base,
                    const int a_dir,
                    const int a_index,
                    const int a_comp) const
    {
      CH_TIME("BSplineInterp::sliceArr");
      ArrSlice slc;
      slc.dir = a_dir;
      slc.size = a_base.size[a_dir];
      slc.stride = a_base.stride[a_dir];
      // go from index to location
      // int slIdx = 0;
      // int idx = a_index;
      // for (int d=0; d!=Dim; d++)
      //   {
      //     if (d != a_dir)
      //       {
      //         slIdx += (idx%a_base.size[d])*a_base.stride[d];
      //         idx = idx/a_base.size[d];
      //       }
      //   }
      // CH_assert(slIdx == m_sIdx[a_dir][a_index]);
      // slIdx += a_comp*product(a_base.size);
      const int slIdx = m_sIdx[a_dir][a_index] +
        a_comp*stc::product(a_base.size);
      slc.data = &a_base.data[slIdx];
      //pout () << slIdx << " " << slc.stride << " " << slc.size << endl;
      return slc;
    }
};

//--"Re" define some constexpr statics

template<int Dim, int Order>
constexpr typename BSplineInterp<Dim, Order>::ISVec
BSplineInterp<Dim, Order>::c_esize;
template<int Dim, int Order>
constexpr typename BSplineInterp<Dim, Order>::ISVec
BSplineInterp<Dim, Order>::c_estride;


#include "NamespaceFooter.H"
#endif
