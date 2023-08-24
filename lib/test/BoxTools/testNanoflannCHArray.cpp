#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <nanoflann.H>

#include "SPACE.H"
#include "REAL.H"
#include "IntVect.H"
#include "RealVect.H"
#include "CHArray.H"
#include "BoxIterator.H"
#include "CH_Timer.H"
#include "parstream.H"
#if CXXSTD>=14
#include "StcVector.H"
#endif

void dump_mem_usage();

// Our array type
using CHarr = CHArray<Real, SpaceDim+1>;


// And this is the "dataset to kd-tree" adaptor class:
struct CHadaptor
{
  const CHarr& m_data; // A const ref to the data set origin

  /// The constructor that sets the data set source
  CHadaptor(const CHarr& a_data)
    :
    m_data(a_data)
    { }

  /// CRTP helper method
  inline const CHarr& data() const
    {
      return m_data;
    }
  
  // Must return the number of data points
  inline size_t kdtree_get_point_count() const
    {
      auto sizeV(m_data.sizeVector());
      sizeV[SpaceDim] = 1; // don't count the number of components
      auto size = stc::product(sizeV);
      CH_assert(size >= 0);
      return size;
    }

  // Returns the dim'th component of the idx'th point in the class:
  inline Real kdtree_get_pt(const size_t idx, const size_t dim) const
    {
      CH_assert(dim >= 0);
      CH_assert(dim < SpaceDim);
      auto leng = kdtree_get_point_count();
      CH_assert(idx < leng);
      return m_data.begin()[idx + dim*leng];
    }

  // method for accessing elements, not required as part of adaptor
  inline RealVect get_pt(const size_t idx) const
    {
      RealVect sol;
      for (auto dir : EachDir)
        {
          sol[dir] = kdtree_get_pt(idx, dir);
        }
      return sol;
    }
  
  // Optional bounding-box computation: return false to determine to a standard bbox computation loop.
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};


void generatePoints(CHarr& data, unsigned type = 0)
{
  pout() << "Generating ";
  if (type == 0)
    {
      pout() << "random points [0, 1] ... ";
      for (auto i = data.begin(); i != data.end(); ++i)
        {
          *i = ((std::rand() % 1000) / Real(1000));
        }
    }
  else if (type == 1)
    {
      pout() << "linear mapping [0, 1] ... ";
      IntVect isize(data.sizeVector());
      Box arrayBox(IntVect::Zero, isize-IntVect::Unit);
      BoxIterator bit(arrayBox);
      for (bit.begin(); bit.ok(); ++bit)
        {
          for (auto dir : EachDir)
            {
              data(bit(), dir) = bit()[dir]/Real(isize[dir]);
            }
        }
    }
  pout() << "done" << std::endl;
  pout() << data << std::endl;
}


void testNanoflann(const size_t N)
{
  CH_TIMERS("testNanoflann");
  
  // Generate data points:
  CHarr data;
  IntVect size(11*IntVect::Unit);
  data.define(size, SpaceDim);
  generatePoints(data, 1);

  // Points to evaluate
  std::vector<RealVect> evalPts {
    RealVect::Zero,
      0.5*RealVect::Unit,
      RealVect{0.4, 0.75},
      RealVect::Unit};

  // Define the adaptor for the point cloud
  const CHadaptor arrayToKD(data);
  pout() << "Array size " << arrayToKD.kdtree_get_point_count() << std::endl;

  // construct a kd-tree index
  dump_mem_usage();
  using CH_KDtree = 
    nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<Real, CHadaptor >,
                                        CHadaptor,
                                        SpaceDim>;
  CH_TIMER("tree construction", t1);
  CH_START(t1);
  CH_KDtree index(SpaceDim,
                  arrayToKD,
                  nanoflann::KDTreeSingleIndexAdaptorParams(50 /* max leaf */) );
  index.buildIndex();
  CH_STOP(t1);
  dump_mem_usage();

  // do a knn search
  CH_TIMER("tree eval", t2);
  CH_START(t2);
  const size_t num_results = 1;
  for (auto pt : evalPts)
    {
      size_t ret_index;
      Real out_dist_sqr;
      nanoflann::KNNResultSet<Real> resultSet(num_results);
      resultSet.init(&ret_index, &out_dist_sqr);
      index.findNeighbors(resultSet, pt.dataPtr(), nanoflann::SearchParams());
      //index.knnSearch(query, indices, dists, num_results, mrpt_flann::SearchParams());
      // extract value
      RealVect sol = arrayToKD.get_pt(ret_index);
      // the index
      Box nodeBox(IntVect::Zero, IntVect(data.sizeVector())-IntVect::Unit);
      BoxIterator bit(nodeBox);
      IntVect KDidx(bit.at(ret_index));
      // print
      pout() << "findNeighbors(nn=" << num_results << "): ";
      pout() << " index = " << ret_index
             << " dist_sqr = " << out_dist_sqr
             << " value = " << pt
             << " neighbor = " << sol
             << " neighbor = " << data(KDidx, 0) << " " << data(KDidx, 1)
             << std::endl;
    }
  CH_STOP(t2);
}

int main()
{
  // Randomize Seed
  std::srand((unsigned int)time(NULL));
  testNanoflann(10E6);
  return 0;
}

void dump_mem_usage()
{
  FILE* f=fopen("statemem","rt");
  if (!f) return;
  char str[300];
  size_t n=fread(str,1,200,f);
  str[n]=0;
  printf("MEM: %s\n",str);
  fclose(f);
}
