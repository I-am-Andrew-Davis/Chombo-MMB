#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstring>
#include <set>
#include <map>
#include <functional>
#include <utility>

#if CXXSTD>=14
#define D_DECL6(a,b,c,d,e,f) a,b,c
#include "StcVector.H"
#endif
#include "CH_Timer.H"
#include "parstream.H"
#include "ComputeCache.H"

#include "UsingBaseNamespace.H"

/// Prototypes:

void
parseTestOptions(int argc, char* argv[]);

int
testCache();

/// Global variables for handling output:
static const char *pgmname = "testCache" ;
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions(argc, argv);

  if (verbose)
    pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl;

  ///
  // Run the tests
  ///
  int ret = testCache();

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << std::endl;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)"
             << std::endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

#if CXXSTD>=14
// Dummy IntVect class with SpaceDim=3
constexpr int c_SpaceDim = 3;
using IntVect = stc::IVec<c_SpaceDim>;
constexpr auto IntVect_zero = stc::make_IVec<c_SpaceDim>::zero();
constexpr auto IntVect_unit = stc::make_IVec<c_SpaceDim>::unit();

namespace CH_Hash
{

/// Remark that IntVect can be byte-hashed and give bytes to hash
template <>
struct isHashable<IntVect> : public std::true_type
{
  static_assert(sizeof(IntVect) == c_SpaceDim*sizeof(int),
                "Cannot guarantee that IntVect is not padded");
  static constexpr int c_hashSize = sizeof(IntVect);
};

}
#endif

int testCache()
{
  CH_TIMERS("testCache");
  int status = 0;

  // Set 1: Use lambda to construct double from IntVect
  {
    stc::RVec<3> results;
    auto cache = make_ComputeCacheUMap<IntVect, double>(
      []
      (const IntVect& a_iv)
      {
        return a_iv.sum() + 0.5;
      });
    results[0] = cache(IntVect{0,1,2});
    results[1] = cache(IntVect{0,1,2});
    results[2] = cache(IntVect{1,2,3});
    if (results != stc::RVec<3>{ 3.5, 3.5, 6.5 }) ++status;
    if (verbose)
      {
        pout() << "Set 1: " << results << std::endl;
      }
  }

  // Set 2: Use std::function containing lambda to construct double from IntVect
  {
    stc::RVec<3> results;
    std::function<double(const IntVect&)> func(
      []
      (const IntVect& a_iv) -> double
      {
        return a_iv.sum() + 0.5;
      });
    auto cache = make_ComputeCacheUMap<IntVect, double>(func);
    results[0] = cache(IntVect{0,1,2});
    results[1] = cache(IntVect{0,1,2});
    results[2] = cache(IntVect{1,2,3});
    if (results != stc::RVec<3>{ 3.5, 3.5, 6.5 }) ++status;
    if (verbose)
      {
        pout() << "Set 2: " << results << std::endl;
      }
    // using X = decltype(cache)::maker_type;
    // pout() << typeid(X).name() << std::endl;
  }

  // Set 3: Same as above but not using maker.  This way you can change the
  // hash, etc.
  {
    stc::RVec<3> results;
    std::function<double(const IntVect&)> func(
      []
      (const IntVect& a_iv) -> double
      {
        return a_iv.sum() + 0.5;
      });
    ComputeCacheUMap<IntVect,
                     double,
                     decltype(make_ComputeCacheUMap<IntVect, double>(func))
                     ::maker_type>
      cache(func);
    results[0] = cache(IntVect{0,1,2});
    results[1] = cache(IntVect{0,1,2});
    results[2] = cache(IntVect{1,2,3});
    if (results != stc::RVec<3>{ 3.5, 3.5, 6.5 }) ++status;
    if (verbose)
      {
        pout() << "Set 3: " << results << std::endl;
      }
  }

  // Set 4: Same as set 2 but using full maker
  {
    stc::RVec<3> results;
    std::function<double(const IntVect&)> func(
      []
      (const IntVect& a_iv) -> double
      {
        return a_iv.sum() + 0.5;
      });
    auto cache = make_ComputeCacheUMap<IntVect,
                                       double,
                                       CH_Hash::google_CityHash<IntVect>,
                                       std::equal_to<IntVect>>
      (func);
    results[0] = cache(IntVect{0,1,2});
    results[1] = cache(IntVect{0,1,2});
    results[2] = cache(IntVect{1,2,3});
    if (results != stc::RVec<3>{ 3.5, 3.5, 6.5 }) ++status;
    if (verbose)
      {
        pout() << "Set 4: " << results << std::endl;
      }
  }

  class Box
  {
  public:
    Box(const IntVect& a_lo, const IntVect& a_hi)
      :
      m_lo(a_lo),
      m_hi(a_hi)
      { }
    double op(const IntVect a_iv)
      {
        return stc::sum(m_lo + m_hi + a_iv) + 0.5;
      }
    IntVect m_lo;
    IntVect m_hi;
  };

  // Set 5: Using lambda to access class member function
  {
    stc::RVec<3> results;
    Box box(IntVect_unit, 2*IntVect_unit);
    auto cache = make_ComputeCacheUMap<IntVect, double>(
      [&]
      (const IntVect& a_iv)
      {
        return box.op(a_iv);
      });
    results[0] = cache(IntVect{0,1,2});
    results[1] = cache(IntVect{0,1,2});
    results[2] = cache(IntVect{1,2,3});
    if (results != stc::RVec<3>{ 12.5, 12.5, 15.5 }) ++status;
    if (verbose)
      {
        pout() << "Set 5: " << results << std::endl;
      }
  }

  // Set 6: Use std::function to access class member function
  // In this case, you must provide the specific arguments for the maker which
  // incudes both the box object and the argument to its member function.
  {
    stc::RVec<3> results;
    Box box(IntVect_unit, 2*IntVect_unit);
    std::function<double(Box&, const IntVect&)> func(&Box::op);
    // The above is an lvalue and would be stored as a reference.  To store it
    // internally in the cache, we can move it as an rvalue.
    auto cache = make_ComputeCacheUMap<IntVect, double>(std::move(func));
    // Try with and without the std::move to func.
    // If the argument is passed as an rvalue, it is stored as a value.
    // If the argument is passed as an lvalue, it is stored as a reference.
    if (verbose)
      {
        using X = decltype(cache)::maker_type;
        std::cout << "Argument: Is rvalue ref: "
                  << std::is_rvalue_reference<X&&>::value << std::endl;
        std::cout << "Argument: Is lvalue ref: "
                  << std::is_lvalue_reference<X&&>::value << std::endl;
        std::cout << "Storage : Is rvalue ref: "
                  << std::is_rvalue_reference<X>::value << std::endl;
        std::cout << "Storage : Is lvalue ref: "
                  << std::is_lvalue_reference<X>::value << std::endl;
      }
    results[0] = cache(IntVect{0,1,2}, box, IntVect{0,1,2});
    results[1] = cache(IntVect{0,1,2}, box, IntVect{0,1,2});
    results[2] = cache(IntVect{1,2,3}, box, IntVect{1,2,3});
    if (results != stc::RVec<3>{ 12.5, 12.5, 15.5 }) ++status;
    if (verbose)
      {
        pout() << "Set 6: " << results << std::endl;
      }
  }

  return status;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions(int argc, char* argv[])
{
  for (int i = 1; i < argc; ++i)
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp( argv[i], "-v", 3 ) == 0)
            {
              verbose = true;
              // argv[i] = "";
            }
          else if (strncmp( argv[i], "-q", 3 ) == 0)
            {
              verbose = false;
              // argv[i] = "";
            }
          else
            {
              break;
            }
        }
    }
  return;
}
