#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <cstring>

#include "parstream.H"
#include "SPMD.H"
#include "DynArray.H"

#ifdef CH_GPU
  #include "testDynArray_CUX.H"
  #include "CudaDriver.H"
#endif

#include "UsingBaseNamespace.H"

/// Prototypes:

Real BaseFabRealSetVal = BASEFAB_REAL_SETVAL;

void
parseTestOptions(int argc, char* argv[]);

int
testDynArray1();

/// Global variables for handling output:
static const char *pgmname = "testDynArray";
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
    {
      pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl;
    }

#ifdef CH_GPU
  CH_Cuda::g_device.initialize(2*(int)verbose);
#endif

  ///
  // Run the tests
  ///
  int ret = testDynArray1();

  if (ret == 0)
    {
      if (verbose)
        {
          pout() << indent << pgmname << " passed all tests" << std::endl;
        }
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)"
             << std::endl ;
    }

#ifdef CH_GPU
  CH_Cuda::g_device.finalize();
#endif
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

// Dummy box class
class Box
{
public:
  Box()
    :
    m_lo(0),
    m_hi(-1)
    { }
  Box(const int a_lo, const int a_hi)
    :
    m_lo(a_lo),
    m_hi(a_hi)
    { }
  int m_lo;
  int m_hi;
};

// These macros ease writing of the nested loops.
#define BOUNDSLOOP3(lb, ub, x)                                          \
  for (int x##2 = (lb)[2],                                              \
         x##2_end = (ub)[2] + 1,                                        \
         x##1_end = (ub)[1] + 1,                                        \
         x##0_end = (ub)[0] + 1;                                        \
       x##2 != x##2_end; ++ x##2)                                       \
    for (int x##1 = (lb)[1]; x##1 != x##1_end; ++ x##1)                 \
      for (int x##0 = (lb)[0]; x##0 != x##0_end; ++ x##0)

#define BOUNDSLOOP2(lb, ub, x)                                          \
  for (int x##1 = (lb)[1],                                              \
         x##1_end = (ub)[1] + 1,                                        \
         x##0_end = (ub)[0] + 1;                                        \
       x##1 != x##1_end; ++ x##1)                                       \
      for (int x##0 = (lb)[0]; x##0 != x##0_end; ++ x##0)

#define SIZELOOP3(arr, x)                                               \
  for (int x##2 = 0,                                                    \
         x##2_end = (arr).size(2),                                      \
         x##1_end = (arr).size(1),                                      \
         x##0_end = (arr).size(1);                                      \
       x##2 != x##2_end; ++ x##2)                                       \
    for (int x##1 = 0; x##1 != x##1_end; ++ x##1)                       \
      for (int x##0 = 0; x##0 != x##0_end; ++ x##0)

int testDynArray1()
{
  int status = 0;

//--Load the module from the GPU (has functions for several tests)

#ifdef CH_GPU
  /// Handles to module global functions for this application
  enum class App_Global
    {
      arrayTestRW = 0,
      arrayTestAliasThread,
      arrayTestAliasShared,
      arrayTestRawThread,
      arrayTestRawShared,
      arrayTestAlloc
    };

  /// Handles to module symbols for this application
  enum class App_Symbol
    {
      g_verbose = 0,
      g_stat,
      g_devicePtr
    };

  // Configure the functions and symbols for this application (test)
  CH_Cuda::g_device.configureChomboLibFunctions(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceFunctionConfig> config;
        config.reserve(6);
        // All kernels use this block size
        dim3 numThread(4, 4, 4);
        config.emplace_back("arrayTestRW", numThread);
        config.emplace_back("arrayTestAliasThread", numThread);
        config.emplace_back("arrayTestAliasShared", numThread);
        config.emplace_back("arrayTestRawThread", numThread);
        config.emplace_back("arrayTestRawShared", numThread);
        config.emplace_back("arrayTestAlloc", numThread);
        return config;
      });

  CH_Cuda::g_device.configureChomboLibSymbols(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceSymbolConfig> config;
        config.reserve(3);
        config.emplace_back("g_verbose");
        config.emplace_back("g_stat");
        config.emplace_back("g_devicePtr");
        return config;
      });

  // Now load the module which automatically loads the functions and symbols
  CH_Cuda::g_device.loadApplicationModule("testDynArray", testDynArray);
#endif

//--Set 1: Test construction (array format)

  {
    int lstat = 0;
    dyn::Array<Real, 3> arrA(4, 3, 2, 0.5);
    for (int k = 0; k != 2; ++k)
      for (int j = 0; j != 3; ++j)
        for (int i = 0; i != 4; ++i)
          {
            if (arrA(i, j, k) != 0.5) ++lstat;
          }
    dyn::Array<Real, 3> arrB(4, 3, 2);
    for (int k = 0; k != 2; ++k)
      for (int j = 0; j != 3; ++j)
        for (int i = 0; i != 4; ++i)
          {
            if (arrB(i, j, k) != BaseFabRealSetVal) ++lstat;
          }
    dyn::Array<Box, 3> arrC(4, 3, 2);
    for (int k = 0; k != 2; ++k)
      for (int j = 0; j != 3; ++j)
        for (int i = 0; i != 4; ++i)
          {
            if (arrC[k][j][i].m_lo != 0 || arrC(i, j, k).m_hi != -1) ++lstat;
          }
    dyn::Array<Box, 3> arrD(4, 3, 2, 10, 12);
    for (int k = 0; k != 2; ++k)
      for (int j = 0; j != 3; ++j)
        for (int i = 0; i != 4; ++i)
          {
            if (arrD[k][j][i].m_lo != 10 || arrD(i, j, k).m_hi != 12) ++lstat;
          }
    stc::Vector<int, 3> lb{ -1, 0, 1 };
    stc::Vector<int, 3> ub{  2, 2, 2 };
    dyn::Array<Real, 3> arrE(lb, ub, 1.4);
    if (arrE.size() != 24) ++lstat;
    if (arrE.size(0) != 4) ++lstat;
    if (arrE.size(1) != 3) ++lstat;
    if (arrE.size(2) != 2) ++lstat;
#ifdef SHAPE_ARRAY_DEBUG
    if (arrE.lowerBound(0) != -1) ++lstat;
    if (arrE.lowerBound(1) != 0) ++lstat;
    if (arrE.lowerBound(2) != 1) ++lstat;
    if (arrE.upperBound(0) != 2) ++lstat;
    if (arrE.upperBound(1) != 2) ++lstat;
    if (arrE.upperBound(2) != 2) ++lstat;
#endif
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrE[i2][i1][i0] != 1.4) ++lstat;
      }
    int c = 0;
    for (const Real& val : arrE)
      {
        ++c;
        if (val != 1.4) ++lstat;
      }
    if (c != arrE.size()) ++lstat;
    if (lstat != 0)
      {
        pout() << "Failure in set 1" << std::endl;
        status += lstat;
      }
  }

//--Set 2: Test move construction and assignment (array format)

  {
    int lstat = 0;
    stc::Vector<int, 3> lb{ -4, -1, 1 };
    stc::Vector<int, 3> ub{ -1,  1, 2 };
    dyn::Array<Box, 3> arrA(lb, ub, 1, 2);
    dyn::Array<Box, 3> arrB(std::move(arrA));
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrB(stc::Vector<int, 3>{ i0, i1, i2 }).m_lo != 1 ||
            arrB(i0, i1, i2).m_hi != 2) ++lstat;
      }
    if (arrA.data() != nullptr) ++lstat;
    dyn::Array<Box, 3> arrC(2, 2, 2, -1, -2);
    arrC = std::move(arrB);
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrC(stc::Vector<int, 3>{ i0, i1, i2 }).m_lo != 1 ||
            arrC(i0, i1, i2).m_hi != 2) ++lstat;
      }
    if (arrB.data() != nullptr) ++lstat;
    if (lstat != 0)
      {
        pout() << "Failure in set 2" << std::endl;
        status += lstat;
      }
  }

//--Set 3: Test construction (alias format)

  {
    int lstat = 0;
    stc::Vector<int, 3> lb{ -4, -1, 1 };
    stc::Vector<int, 3> ub{ -1,  1, 2 };
    dyn::Array<Box, 3> arrA(lb, ub, 1, 2);
    dyn::Array<Box, 3> arrB(dyn::alias_format{}, arrA.data(),
                            arrA.size(0), arrA.size(1), arrA.size(2));
    SIZELOOP3(arrB, i)
      {
        if (arrB(stc::Vector<int, 3>{ i0, i1, i2 }).m_lo != 1 ||
            arrB(i0, i1, i2).m_hi != 2) ++lstat;
      }
    // Test alias with bounds
    dyn::Array<Box, 3> arrC(dyn::alias_format{}, arrA.data(), lb, ub);
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrC(stc::Vector<int, 3>{ i0, i1, i2 }).m_lo != 1 ||
            arrC(i0, i1, i2).m_hi != 2) ++lstat;
      }
    // Test alias with bounds
    dyn::Array<Box, 3> arrD(dyn::alias_format{}, arrA);
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrD(stc::Vector<int, 3>{ i0, i1, i2 }).m_lo != 1 ||
            arrD(i0, i1, i2).m_hi != 2) ++lstat;
      }
    // Test some modifications
    arrD(-3, 0, 2).m_lo = 3;
    arrD(-2, 0, 1).m_hi = 4;
    if (arrA(-3, 0, 2).m_lo != 3) ++lstat;
    if (arrA(-3, 0, 2).m_hi != 2) ++lstat;
    if (arrA(-2, 0, 1).m_lo != 1) ++lstat;
    if (arrA(-2, 0, 1).m_hi != 4) ++lstat;
    if (arrB(1, 1, 1).m_lo != 3) ++lstat;  // arrB is zero-based
    if (arrB(1, 1, 1).m_hi != 2) ++lstat;
    if (arrB(2, 1, 0).m_lo != 1) ++lstat;
    if (arrB(2, 1, 0).m_hi != 4) ++lstat;
    if (arrC(-3, 0, 2).m_lo != 3) ++lstat;
    if (arrC(-3, 0, 2).m_hi != 2) ++lstat;
    if (arrC(-2, 0, 1).m_lo != 1) ++lstat;
    if (arrC(-2, 0, 1).m_hi != 4) ++lstat;
    if (lstat != 0)
      {
        pout() << "Failure in set 3" << std::endl;
        status += lstat;
      }
  }

//--Set 4: Test construction (raw format)

  {
    int lstat = 0;
    stc::Vector<int, 3> lb{ -4, -1, 1 };
    stc::Vector<int, 3> ub{ -1,  1, 2 };
    dyn::Array<Box, 3> arrA(lb, ub, 1, 2);
    // Note: you wouldn't want to do this if Box was not trivial because the
    // destructor will be called multiple times.  Raw should not be used with
    // constructed memory.
    dyn::Array<Box, 3> arrB(dyn::raw_format{}, arrA.data(),
                            arrA.size(0), arrA.size(1), arrA.size(2),
                            -1, -2);
    SIZELOOP3(arrB, i)
      {
        if (arrB(stc::Vector<int, 3>{ i0, i1, i2 }).m_lo != -1 ||
            arrB(i0, i1, i2).m_hi != -2) ++lstat;
      }
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrA(stc::Vector<int, 3>{ i0, i1, i2 }).m_lo != -1 ||
            arrA(i0, i1, i2).m_hi != -2) ++lstat;
      }
    // Test raw with bounds
    dyn::Array<Box, 3> arrC(dyn::raw_format{}, arrA.data(), lb, ub,
                            3, 4);
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrC(stc::Vector<int, 3>{ i0, i1, i2 }).m_lo != 3 ||
            arrC(i0, i1, i2).m_hi != 4) ++lstat;
      }
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrA(stc::Vector<int, 3>{ i0, i1, i2 }).m_lo != 3 ||
            arrA(i0, i1, i2).m_hi != 4) ++lstat;
      }
    if (lstat != 0)
      {
        pout() << "Failure in set 4" << std::endl;
        status += lstat;
      }
  }

//--Set 5: Test define (array format)

  {
    int lstat = 0;
    dyn::Array<Box, 3> arrA(stc::Vector<int, 3>{ -4, -1, 1 },
                            stc::Vector<int, 3>{ -1,  1, 2 }, 1, 2);
    arrA.define(3, 3, 3, -1, -2);
    SIZELOOP3(arrA, i)
      {
        if (arrA[i2][i1][i0].m_lo != -1 ||
            arrA(i0, i1, i2).m_hi != -2) ++lstat;
      }
    arrA.define(2, 2, 2);
    SIZELOOP3(arrA, i)
      {
        if (arrA[i2][i1][i0].m_lo != 0 ||
            arrA(i0, i1, i2).m_hi != -1) ++lstat;
      }
    stc::Vector<int, 3> lb{ -1, -1, 2 };
    stc::Vector<int, 3> ub{  2,  2, 5 };
    arrA.define(lb, ub, 3, 4);
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrA[i2][i1][i0].m_lo != 3 ||
            arrA(i0, i1, i2).m_hi != 4) ++lstat;
      }
    if (lstat != 0)
      {
        pout() << "Failure in set 5" << std::endl;
        status += lstat;
      }
  }

//--Set 6: Test define (alias format)

  {
    int lstat = 0;
    dyn::Array<Box, 3> arrA(3, 3, 3, 1, 2);
    for (int j = 0; j != 3; ++j)
      for (int i = 0; i != 3; ++i)
        {
          arrA(i, j, 1).m_lo = -1;
          arrA(i, j, 1).m_hi = -2;
        }
    dyn::Array<Box, 3> arrB;
    arrB.defineAlias(&arrA(0, 0, 1), 3, 3, 1);
    SIZELOOP3(arrB, i)
      {
        if (arrB[i2][i1][i0].m_lo != -1 ||
            arrB(i0, i1, i2).m_hi != -2) ++lstat;
        arrB(i0, i1, i2).m_lo = -3;
      }
    for (int j = 0; j != 3; ++j)
      for (int i = 0; i != 3; ++i)
        {
          if (arrA(i, j, 1).m_lo != -3) ++lstat;
        }
    // Alias with resize
    stc::Vector<int, 2> lb{ -1, 0 };
    stc::Vector<int, 2> ub{  1, 2 };
    dyn::Array<Box, 2> arrC(1, 1, 1, 2);
    arrC.defineAlias(&arrA(0, 0, 1), lb, ub);
    BOUNDSLOOP2(lb, ub, i)
      {
        if (arrC[i1][i0].m_lo != -3 ||
            arrC(i0, i1).m_hi != -2) ++lstat;
      }
    // Alias full array
    dyn::Array<Box, 2> arrD(1, 1, 1, 2);
    arrD.defineAlias(arrC);
    BOUNDSLOOP2(lb, ub, i)
      {
        if (arrD[i1][i0].m_lo != -3 ||
            arrD(i0, i1).m_hi != -2) ++lstat;
        arrD(i0, i1).m_lo = -4;
      }
    for (int j = 0; j != 3; ++j)
      for (int i = 0; i != 3; ++i)
        {
          if (arrA(i, j, 1).m_lo != -4) ++lstat;
        }
    if (lstat != 0)
      {
        pout() << "Failure in set 6" << std::endl;
        status += lstat;
      }
  }

//--Set 7: Test define (raw format)

  {
    int lstat = 0;
    unsigned char* data = new unsigned char[3*3*3*sizeof(Box)];
    dyn::Array<Box, 3> arrA(1, 1, 1, 3, 4);
    arrA.defineRaw(data, 3, 3, 3, 1, 2);
    SIZELOOP3(arrA, i)
      {
        if (arrA[i2][i1][i0].m_lo != 1 ||
            arrA(i0, i1, i2).m_hi != 2) ++lstat;
      }
    // Raw with resize
    stc::Vector<int, 2> lb{ -1, 0 };
    stc::Vector<int, 2> ub{  1, 2 };
    dyn::Array<Box, 2> arrB(1, 1, 3, 4);
    arrB.defineRaw(data + 3*3*sizeof(Box), lb, ub, -1, -2);
    BOUNDSLOOP2(lb, ub, i)
      {
        if (arrB[i1][i0].m_lo != -1 ||
            arrB(i0, i1).m_hi != -2) ++lstat;
      }
    // Still constructs using default constructor for box
    dyn::Array<Box, 2> arrC(1, 1, 3, 4);
    arrC.defineRaw(data + 3*3*sizeof(Box), lb, ub);
    BOUNDSLOOP2(lb, ub, i)
      {
        if (arrC[i1][i0].m_lo !=  0 ||
            arrC(i0, i1).m_hi != -1) ++lstat;
      }
    // Test the original
    for (int j = 0; j != 3; ++j)
      for (int i = 0; i != 3; ++i)
        {
          if (arrA[0][j][i].m_lo !=  1 ||
              arrA(i, j, 0).m_hi !=  2) ++lstat;
          if (arrA[1][j][i].m_lo !=  0 ||
              arrA(i, j, 1).m_hi != -1) ++lstat;
          if (arrA[2][j][i].m_lo !=  1 ||
              arrA(i, j, 2).m_hi !=  2) ++lstat;
        }
    delete[] data;
    if (lstat != 0)
      {
        pout() << "Failure in set 7" << std::endl;
        status += lstat;
      }
  }

//--Set 8: Test initialization

  {
    int lstat = 0;
    stc::Vector<int, 3> lb{ -4, -1, 0 };
    stc::Vector<int, 3> ub{ -1,  1, 2 };
    dyn::Array<Box, 3> arrA(lb, ub, 1, 2);
    arrA = Box(3, 4);
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrA(i0, i1, i2).m_lo != 3 ||
            arrA(i0, i1, i2).m_hi != 4) ++lstat;
      }
    dyn::Array<Box, 3> arrB(dyn::alias_format{},
                            &arrA(-4, -1, 1),
                            arrA.size(0), arrA.size(1), 1);
    arrB = Box(5, 6);
    BOUNDSLOOP2(lb, ub, i)
      {
        if (arrA(i0, i1, 0).m_lo != 3 || arrA(i0, i1, 0).m_hi != 4) ++lstat;
        if (arrA(i0, i1, 1).m_lo != 5 || arrA(i0, i1, 1).m_hi != 6) ++lstat;
        if (arrA(i0, i1, 2).m_lo != 3 || arrA(i0, i1, 2).m_hi != 4) ++lstat;
      }
    if (lstat != 0)
      {
        pout() << "Failure in set 8" << std::endl;
        status += lstat;
      }
  }

//--Set 9: Most [] and () indexing operators were tested in the above.  Test
//--indexing with vector plus additional indices.

  {
    int lstat = 0;
    stc::Vector<int, 4> lb{ -1, -1, 0, 0 };
    stc::Vector<int, 4> ub{  2,  2, 2, 1 };
    dyn::Array<Box, 4> arrA(lb, ub);
    for (int i3 = lb[3]; i3 <= lb[3]; ++i3)
      BOUNDSLOOP3(lb, ub, i)
        {
          arrA[i3][i2][i1][i0].m_lo = i0 + i1 + i2 + i3;
          arrA[i3][i2][i1][i0].m_hi = -(i0 + i1 + i2 + i3);
        }
    for (int i3 = lb[3]; i3 <= lb[3]; ++i3)
      for (int i2 = lb[2]; i2 <= lb[2]; ++i2)
        {
          stc::Vector<int, 2> iv(i2, i3);
          for (int i1 = lb[1]; i1 <= lb[1]; ++i1)
            for (int i0 = lb[0]; i0 <= lb[0]; ++i0)
              {
                if (arrA(i0, i1, iv).m_lo != i0 + i1 + stc::sum(iv)) ++lstat;
                if (arrA(i0, i1, iv).m_hi != -(i0 + i1 + stc::sum(iv))) ++lstat;
              }
        }
    for (int i3 = lb[3]; i3 <= lb[3]; ++i3)
      BOUNDSLOOP3(lb, ub, i)
        {
          stc::Vector<int, 2> iv(i0, i1);
          if (arrA(iv, i2, i3).m_lo != stc::sum(iv) + i2 + i3) ++lstat;
          if (arrA(iv, i2, i3).m_hi != -(stc::sum(iv) + i2 + i3)) ++lstat;
        }
    if (lstat != 0)
      {
        pout() << "Failure in set 9" << std::endl;
        status += lstat;
      }
  }

//--Set 10: Test shift

  {
    int lstat = 0;
    dyn::Array<Box, 3> arrA(3, 3, 3, -1, -2);
    arrA.shift(-1, 0, 1);
    stc::Vector<int, 3> lb{ -1, 0, 1 };
    stc::Vector<int, 3> ub{  1, 2, 3 };
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrA(i0, i1, i2).m_lo != -1) ++lstat;
        if (arrA(i0, i1, i2).m_hi != -2) ++lstat;
      }
    if (&arrA(lb) != arrA.data()) ++lstat;
    if (&arrA(lb) != arrA.begin()) ++lstat;
    // Shift again with a vector
    arrA.shift(lb);
    ub += lb;
    lb += lb;
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrA[i2][i1][i0].m_lo != -1) ++lstat;
        if (arrA(i0, i1, i2).m_hi != -2) ++lstat;
      }
    if (&arrA(lb) != arrA.data()) ++lstat;
    if (&arrA(lb) != arrA.begin()) ++lstat;

    // Repeat with a row-ordered alias
    dyn::Array<Box, 3, int, shape::row_ordered> arrB(
      dyn::alias_format{}, arrA.data(), 3, 3, 3);
    arrB.shift(1, 0, -1);
    lb = stc::Vector<int, 3>{ 1, 0, -1 };
    ub = stc::Vector<int, 3>{ 3, 2,  1 };
    BOUNDSLOOP3(lb, ub, i)
      {
        // Note: the indices in bounds loop are backwards (i2 is really i0), so
        // the order appears the same here using (), even though arrB is
        // row-ordered.  See differences in operator[] from above.
        if (arrB[i0][i1][i2].m_lo != -1) ++lstat;
        if (arrB(i0, i1, i2).m_hi != -2) ++lstat;
      }
    if (&arrB(lb) != arrB.data()) ++lstat;
    if (&arrB(lb) != arrB.begin()) ++lstat;
    // Shift again with a vector
    arrB.shift(lb);
    ub += lb;
    lb += lb;
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrB[i0][i1][i2].m_lo != -1) ++lstat;
        if (arrB(i0, i1, i2).m_hi != -2) ++lstat;
      }
    if (&arrB(lb) != arrB.data()) ++lstat;
    if (&arrB(lb) != arrB.begin()) ++lstat;
    BOUNDSLOOP3(lb, ub, i)
      {
        if (&arrB[i0][i1][i2] != &arrA[i0][i1][i2]) ++lstat;
        if (&arrB(i0, i1, i2) != &arrA(i2, i1, i0)) ++lstat;
      }
    if (lstat != 0)
      {
        pout() << "Failure in set 10" << std::endl;
        status += lstat;
      }
  }

//--Set 11: Test row ordering

  {
    int lstat = 0;
    stc::Vector<int, 3> lb{ -4, -1, 0 };
    stc::Vector<int, 3> ub{ -1,  1, 2 };
    dyn::Array<Box, 3> arrA(lb, ub);
    int c = 0;
    BOUNDSLOOP3(lb, ub, i)
      {
        arrA[i2][i1][i0].m_lo = c++;
        arrA[i2][i1][i0].m_hi = -(c++);
      }
    dyn::Array<Box, 3, int, shape::row_ordered> arrB(
      dyn::alias_format{}, arrA.data(),
      arrA.size(2), arrA.size(1), arrA.size(0));
    arrB.shift(0, -1, -4);
    BOUNDSLOOP3(lb, ub, i)
      {
        // Bracket indexing always as row ordered
        if (arrA[i2][i1][i0].m_lo != arrB[i2][i1][i0].m_lo) ++lstat;
        if (arrA[i2][i1][i0].m_hi != arrB[i2][i1][i0].m_hi) ++lstat;
        // Parenthesis indexing varies according to order
        if (arrA(i0, i1, i2).m_lo != arrB(i2, i1, i0).m_lo) ++lstat;
        if (arrA(i0, i1, i2).m_hi != arrB(i2, i1, i0).m_hi) ++lstat;
        if (&arrA(i0, i1, i2).m_lo != &arrB(i2, i1, i0).m_lo) ++lstat;
        if (&arrA(i0, i1, i2).m_hi != &arrB(i2, i1, i0).m_hi) ++lstat;
        // Vectors just replace indices
        if (arrA(stc::Vector<int, 2>(i0, i1), i2).m_lo !=
            arrB(i2, stc::Vector<int, 2>(i1, i0)).m_lo) ++lstat;
        if (arrA(i0, i1, stc::Vector<int, 1>(i2)).m_hi !=
            arrB(stc::Vector<int, 3>(i2, i1, i0)).m_hi) ++lstat;
      }
    if (lstat != 0)
      {
        pout() << "Failure in set 11" << std::endl;
        status += lstat;
      }
  }

//--Set 12: These sizes should remain true

  {
    int lstat = 0;
    // The base storage is the pointer (8 bytes) and the offset parameter
    // (8 bytes)
#ifdef SHAPE_ARRAY_DEBUG
    // Debug requires storing lower and upper bounds (8 bytes) per dimension for
    // range checking
#ifdef CH_GPU
    // Saving pointers in the allocator costs an extra 8 bytes
    if (sizeof(dyn::Array<Real, 1>) != 32) ++lstat;
    if (sizeof(dyn::Array<Real, 2>) != 40) ++lstat;
    if (sizeof(dyn::Array<Real, 3>) != 48) ++lstat;
    if (sizeof(dyn::Array<Real, 4>) != 56) ++lstat;
#else
    if (sizeof(dyn::Array<Real, 1>) != 24) ++lstat;
    if (sizeof(dyn::Array<Real, 2>) != 32) ++lstat;
    if (sizeof(dyn::Array<Real, 3>) != 40) ++lstat;
    if (sizeof(dyn::Array<Real, 4>) != 48) ++lstat;
#endif
#else
    // Without range checking, only a stride (4 bytes) is stored per dimension
#ifdef CH_GPU
    // Saving pointers in the allocator costs an extra 8 bytes
    if (sizeof(dyn::Array<Real, 1>) != 32) ++lstat;
    if (sizeof(dyn::Array<Real, 2>) != 32) ++lstat;
    if (sizeof(dyn::Array<Real, 3>) != 40) ++lstat;
    if (sizeof(dyn::Array<Real, 4>) != 40) ++lstat;
#else
    if (sizeof(dyn::Array<Real, 1>) != 24) ++lstat;
    if (sizeof(dyn::Array<Real, 2>) != 24) ++lstat;
    if (sizeof(dyn::Array<Real, 3>) != 32) ++lstat;
    if (sizeof(dyn::Array<Real, 4>) != 32) ++lstat;
#endif
#endif
    if (verbose)
      {
        std::cout << "sizeof(dyn::Array<Real, 1>): "
                  << sizeof(dyn::Array<Real, 1>) << std::endl;
        std::cout << "sizeof(dyn::Array<Real, 2>): "
                  << sizeof(dyn::Array<Real, 2>) << std::endl;
        std::cout << "sizeof(dyn::Array<Real, 3>): "
                  << sizeof(dyn::Array<Real, 3>) << std::endl;
        std::cout << "sizeof(dyn::Array<Real, 4>): "
                  << sizeof(dyn::Array<Real, 4>) << std::endl;
      }
    if (lstat != 0)
      {
        pout() << "Failure in set 12" << std::endl;
        status += lstat;
      }
  }

//--Set 13: Test basic reading/writing on the GPU

#ifdef CH_GPU
  {
    int lstat = 0;

    stc::Vector<int, 3> lb{ -1, 0, 1 };
    stc::Vector<int, 3> ub{  2, 3, 4 };
    dyn::Array<Real, 3> arrA(lb, ub);
    BOUNDSLOOP3(lb, ub, i)
      {
        arrA(i0, i1, i2) = i0 + i1 + i2;
      }

    // Pointer locations
    if (verbose)
      {
        pout() << "Array pointer locations:\n";
        static constexpr const char* ptrAttrLbl[] =
          { "unknown", "host  ", "device" };
        static_assert(CU_MEMORYTYPE_HOST == 1u, "Fix attribute mapping");
        static_assert(CU_MEMORYTYPE_DEVICE == 2u, "Fix attribute mapping");
        unsigned ptrAttr;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            arrA.data())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "arrA on " << ptrAttrLbl[ptrAttr] << "         : " <<
          arrA.data() << std::endl;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            arrA.devicePtr())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "arrA on " << ptrAttrLbl[ptrAttr] << "         : " <<
          arrA.devicePtr() << std::endl;
      }
    
    // Transfer variables to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_verbose,
                                              verbose);
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);
    uintptr_t l_devicePtr = reinterpret_cast<uintptr_t>(arrA.devicePtr());
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_devicePtr,
                                              l_devicePtr);

    arrA.copyToDeviceAsync();
    // Launch the kernel
    dim3 numBlock;
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestRW,
                             numBlock,
                             CH_Cuda::c_defaultStream,
                             arrA);
    arrA.copyToHostAsync();

    // Transfer status to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));
    
    BOUNDSLOOP3(lb, ub, i)
      {
        if (arrA(i0, i1, i2) != i0 + i1 + i2 + 0.5) ++lstat;
      }

    if (lstat != 0)
      {
        pout() << "Failure in set 13" << std::endl;
        status += lstat;
      }
  }
#endif

//--Set 14: Test aliases on the GPU

#ifdef CH_GPU
  {
    int lstat = 0;

    stc::Vector<int, 3> lb{ -1, 0, 1 };
    stc::Vector<int, 3> ub{  2, 3, 4 };
    dyn::Array<Real, 3> arrA(lb, ub, 1.5);

//--Test using dyn::Array per thread

    // Transfer status to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);

    arrA.copyToDeviceAsync();
    // Launch the kernel
    dim3 numBlock;
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestAliasThread,
                             numBlock,
                             CH_Cuda::c_defaultStream,
                             lb, arrA);
    arrA.copyToHostAsync();

    // Transfer status to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));

    for (int k = lb[2]; k <= ub[2]; ++k)
      {
        Real val = 1.5;
        if (k == 0 + lb[2])
          {
            val = 2.5;
          }
        else if (k == 1 + lb[2])
          {
            val = 3.5;
          }
        else if (k == 2 + lb[2])
          {
            val = 4.5;
          }
        for (int j = lb[1]; j <= ub[1]; ++j)
          for (int i = lb[0]; i <= ub[0]; ++i)
            {
              if (arrA(i, j, k) != val) ++lstat;
            }
      }

//--Test using dyn::Array in shared memory

    arrA = 1.5;

    // Transfer status to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);

    arrA.copyToDeviceAsync();
    // Launch the kernel
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestAliasShared,
                             numBlock,
                             CH_Cuda::c_defaultStream,
                             lb, arrA);
    arrA.copyToHostAsync();

    // Transfer status to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));

    for (int k = lb[2]; k <= ub[2]; ++k)
      {
        Real val = 1.5;
        if (k == 0 + lb[2])
          {
            val = 2.5;
          }
        else if (k == 1 + lb[2])
          {
            val = 3.5;
          }
        else if (k == 2 + lb[2])
          {
            val = 4.5;
          }
        for (int j = lb[1]; j <= ub[1]; ++j)
          for (int i = lb[0]; i <= ub[0]; ++i)
            {
              if (arrA(i, j, k) != val) ++lstat;
            }
      }

    if (lstat != 0)
      {
        pout() << "Failure in set 14" << std::endl;
        status += lstat;
      }
  }
#endif

//--Set 15: Test using raw memory on the GPU

#ifdef CH_GPU
  {
    int lstat = 0;
    // Transfer status to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);

    // Launch the kernels
    dim3 numBlock;
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestRawThread,
                             numBlock,
                             CH_Cuda::c_defaultStream);
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestRawShared,
                             numBlock,
                             CH_Cuda::c_defaultStream);

    // Transfer status to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));

    if (lstat != 0)
      {
        pout() << "Failure in set 15" << std::endl;
        status += lstat;
      }
  }
#endif

//--Set 16: Allocations on the GPU

#ifdef CH_GPU
  {
    int lstat = 0;
    // Transfer status to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);

    // Launch the kernel
    dim3 numBlock;
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestAlloc,
                             numBlock,
                             CH_Cuda::c_defaultStream);

    // Transfer status to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));

    if (lstat != 0)
      {
        pout() << "Failure in set 16" << std::endl;
        status += lstat;
      }
  }
#endif

  return status;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp( argv[i] ,"-v" ,3 ) == 0)
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if (strncmp( argv[i] ,"-q" ,3 ) == 0)
            {
              verbose = false ;
              // argv[i] = "" ;
            }
          else
            {
              break ;
            }
        }
    }
  return ;
}
