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

#include "parstream.H"
#include "SPMD.H"

#if CXXSTD>=14
#define D_DECL6(a,b,c,d,e,f) a,b,c
#include "StcVector.H"
#endif
#include "CH_Timer.H"
#include "parstream.H"

void dump_mem_usage();

// This is an exampleof a custom data set class
template <typename T>
struct PointCloud
{
  typedef T coord_t; //!< The type of each coordinate

  struct Point
  {
    T  x,y,z;
  };

  std::vector<Point>  pts;
}; // end of PointCloud

// And this is the "dataset to kd-tree" adaptor class:
template <typename Derived>
struct PointCloudAdaptor
{
  typedef typename Derived::coord_t coord_t;

  const Derived &obj; //!< A const ref to the data set origin

  /// The constructor that sets the data set source
  PointCloudAdaptor(const Derived &obj_) : obj(obj_) { }

  /// CRTP helper method
  inline const Derived& derived() const { return obj; }

  // Must return the number of data points
  inline size_t kdtree_get_point_count() const { return derived().pts.size(); }

  // Returns the dim'th component of the idx'th point in the class:
  // Since this is inlined and the "dim" argument is typically an immediate value, the
  //  "if/else's" are actually solved at compile time.
  inline coord_t kdtree_get_pt(const size_t idx, const size_t dim) const
    {
      if (dim == 0) return derived().pts[idx].x;
      else if (dim == 1) return derived().pts[idx].y;
      else return derived().pts[idx].z;
    }

  // Optional bounding-box computation: return false to default to a standard bbox computation loop.
  //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
  //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

}; // end of PointCloudAdaptor


template <typename T>
void generateRandomPointCloud(PointCloud<T> &point, const size_t N, const T max_range = 10)
{
  std::cout << "Generating "<< N << " point cloud...";
  point.pts.resize(N);
  // for (size_t i = 0; i < N;i++)
  //   {
  //     point.pts[i].x = max_range * (rand() % 1000) / T(1000);
  //     point.pts[i].y = max_range * (rand() % 1000) / T(1000);
  //     point.pts[i].z = max_range * (rand() % 1000) / T(1000);
  //   }
  for (size_t i = 0; i < N;i++)
    {
      point.pts[i].x = max_range * (i/T(N));
      point.pts[i].y = max_range * (i/T(N));
      point.pts[i].z = max_range * (i/T(N));
    }

  std::cout << "done\n";
}

template <typename num_t>
void testNanoflann(const size_t N)
{
  CH_TIMERS("testNanoflann");
  const size_t Dim = 3;

  // Generate data points:
  PointCloud<num_t> cloud;
  generateRandomPointCloud(cloud, N);

  // Points to evaluate
  const size_t n_eval = 7;
  PointCloud<num_t> query_pts;
  generateRandomPointCloud(query_pts, n_eval);

  // Define the adaptor for the point cloud
  typedef PointCloudAdaptor<PointCloud<num_t> > PC2KD;
  const PC2KD pc2kd(cloud); // The adaptor

  // construct a kd-tree index
  dump_mem_usage();
  typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<num_t, PC2KD > ,
                                              PC2KD,
                                              Dim> my_kd_tree_t;
  CH_TIMER("tree construction", t1);
  CH_START(t1);
  my_kd_tree_t index(Dim,
                     pc2kd,
                     nanoflann::KDTreeSingleIndexAdaptorParams(50 /* max leaf */) );
  index.buildIndex();
  CH_STOP(t1);
  dump_mem_usage();

  // do a knn search
  CH_TIMER("tree evals", t2);
  CH_START(t2);
  const size_t num_results = 1;
  for (size_t eval = 0; eval!=n_eval; eval++)
    {
      size_t ret_index;
      num_t out_dist_sqr;
      nanoflann::KNNResultSet<num_t> resultSet(num_results);
      resultSet.init(&ret_index, &out_dist_sqr);
      index.findNeighbors(resultSet, &query_pts.pts[eval].x, nanoflann::SearchParams());
      //index.knnSearch(&query_pts.pts[eval].x, indices, dists, num_results, mrpt_flann::SearchParams());
      pout() << "findNeighbors(nn=" << num_results << "): \t";
      pout() << " index = " << ret_index
             << " dist_sqr = " << out_dist_sqr
             << " value = " << -0 << std::endl;
    }
  CH_STOP(t2);
}

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Randomize Seed
  srand((unsigned int)time(NULL));
  testNanoflann<double>(10E6);
#ifdef CH_MPI
  MPI_Finalize();
#endif
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
