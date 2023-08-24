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
using std::endl;
#include "parstream.H"
#include "MDBitSet.H"
#include "IntVect.H"
#ifdef CH_MPI
#include <mpi.h>
#endif

#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testMDBitSet();

/// Global variables for handling output:
static const char *pgmname = "testMDBitSet" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  if (verbose)
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testMDBitSet() ;

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)" << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


/*******************************************************************************
 *
 * Routine testCHArray
 *
 ******************************************************************************/

int
testMDBitSet()
{
  int status = 0;

//--Experimentation

  using BS = MDBitSet<unsigned, int, unsigned, 8, 8, 16, 16>;
  BS bs;
  // bs[0] = 1;
  // std::cout << bs.count() << std::endl;
  std::cout << bs.rank << std::endl;
  std::cout << bs.size() << std::endl;
  std::cout << bs.NW << std::endl;
  std::cout << bs.dimensions() << std::endl;
  std::cout << "tallyb     : " << bs.c_tallyb << std::endl;
  std::cout << "dimsW      : " << bs.c_dimsW << std::endl;
  std::cout << "dimsWD0inW : " << bs.c_dimsWD0inW << std::endl;
  // std::cout << "strideWallb: " << bs.c_strideWallb << std::endl;
  // std::cout << "strideW    : " << bs.c_strideW << std::endl;
  // std::cout << "strideWin0 : " << bs.c_strideWin0 << std::endl;
  // Need assert that verifies this
  std::cout << stc::max((BS::uix_type)1, bs.c_dimsW).product() << ' ' << bs.NW
            << std::endl;
  BS::IIxVec uv{0, 4, 0};
  std::cout << bs.getb(uv) << std::endl;
  bs.setb(BS::IIxVec{0, 4, 0});
  bs.setb(BS::IIxVec{1, 4, 0});
  bs.setb(BS::IIxVec{2, 4, 0});
  bs.setb(BS::IIxVec{7, 4, 0});
  bs.setb(BS::IIxVec{0, 4, 1});
  std::cout << bs.getb(uv) << std::endl;
  // std::cout << bs.getb(BS::IIxVec{-1, 3, -1}) << std::endl;
  std::cout << bs.getb(BS::UIxVec{7, 3, 0}) << std::endl;
  std::cout << bs.getb(BS::UIxVec{1, 4, 0}) << std::endl;
  std::cout << "count: " << bs.count() << std::endl;
  bs.forEachTrueBit(BS::IIxVec(0),
                    []
                    (const BS::IIxVec& a_iv)
                      {
                        std::cout << a_iv << std::endl;
                      });
  std::cout << "----\n";
  bs.forEachTrueBit(BS::IIxVec(0),
                    BS::IIxVec{1, 0, 0}, BS::IIxVec{2, 15, 15},
                    []
                    (const BS::IIxVec& a_iv)
                      {
                        std::cout << a_iv << std::endl;
                      });
  std::cout << "----\n";
  bs.setAll();
  std::cout << bs.isEmpty() << ' ' << bs.isFull() << std::endl;
  std::cout << "count: " << bs.count() << std::endl;
  bs.clearAll();
  std::cout << bs.isEmpty() << ' ' << bs.isFull() << std::endl;
  std::cout << "count: " << bs.count() << std::endl;

//--Set 1 (constants and word sizes)

  {
    // Check Nb WszB Wszb NW Nbrem
    // unsigned (size 32 D0)
    using BSu32_32 = MDBitSet<unsigned, int, unsigned, 32, 32, 16, 16>;
    static_assert(BSu32_32::rank == 3, "Set1 (u32_32): Invalid rank");
    static_assert(BSu32_32::Nb == 8192, "Set1 (u32_32): Invalid num bits");
    static_assert(BSu32_32::WszB == 4,
                  "Set1 (u32_32): Invalid word size (B)");
    static_assert(BSu32_32::Wszb == 32,
                  "Set1 (u32_32): Invalid word size (b)");
    static_assert(BSu32_32::NW == 256, "Set1 (u32_32): Invalid num words");
    static_assert(BSu32_32::Nbrem == 0,
                  "Set1 (u32_32): Invalid remainder bits");
    static_assert(BSu32_32::Nbsrp == 0,
                  "Set1 (u32_32): Invalid surplus bits");
    static_assert(BSu32_32::c_dimsb == stc::Vector<unsigned, 3>{32, 16, 16},
                  "Set1 (u32_32): Invalid dimensions (u)");
    static_assert(BSu32_32::c_iix_dimsb == stc::Vector<int, 3>{32, 16, 16},
                  "Set1 (u32_32): Invalid dimensions (i)");
    static_assert(BSu32_32::c_strideb == stc::Vector<unsigned, 3>{1, 32, 32*16},
                  "Set1 (u32_32): Invalid bit stride (u)");
    static_assert(BSu32_32::c_iix_strideb == stc::Vector<int, 3>{1, 32, 32*16},
                  "Set1 (u32_32): Invalid bit stride (i)");
    static_assert(BSu32_32::c_D0inWszb == 32,
                  "Set1 (u32_32): Invalid D0inW bit size");
    static_assert(BSu32_32::c_stepD0inWb == stc::Vector<unsigned, 3>{32, 1, 1},
                  "Set1 (u32_32): Invalid D0inW bit steps (u)");
    static_assert(BSu32_32::c_iix_stepD0inWb == stc::Vector<int, 3>{32, 1, 1},
                  "Set1 (u32_32): Invalid D0inW bit steps (i)");
    static_assert(BSu32_32::c_maskD0inW == 0xFFFFFFFF,
                  "Set1 (u32_32): Invalid D0inW mask");
    // unsigned (size 16 D0)
    using BSu32_16 = MDBitSet<unsigned, int, unsigned, 16, 16, 16, 16>;
    static_assert(BSu32_16::rank == 3, "Set1 (u32_16): Invalid rank");
    static_assert(BSu32_16::Nb == 4096, "Set1 (u32_16): Invalid num bits");
    static_assert(BSu32_16::WszB == 4,
                  "Set1 (u32_16): Invalid word size (B)");
    static_assert(BSu32_16::Wszb == 32,
                  "Set1 (u32_16): Invalid word size (b)");
    static_assert(BSu32_16::NW == 128, "Set1 (u32_16): Invalid num words");
    static_assert(BSu32_16::Nbrem == 0,
                  "Set1 (u32_16): Invalid remainder bits");
    static_assert(BSu32_16::Nbsrp == 0,
                  "Set1 (u32_16): Invalid surplus bits");
    static_assert(BSu32_16::c_dimsb == stc::Vector<unsigned, 3>{16, 16, 16},
                  "Set1 (u32_16): Invalid dimensions (u)");
    static_assert(BSu32_16::c_iix_dimsb == stc::Vector<int, 3>{16, 16, 16},
                  "Set1 (u32_16): Invalid dimensions (i)");
    static_assert(BSu32_16::c_strideb == stc::Vector<unsigned, 3>{1, 16, 16*16},
                  "Set1 (u32_16): Invalid bit stride (u)");
    static_assert(BSu32_16::c_iix_strideb == stc::Vector<int, 3>{1, 16, 16*16},
                  "Set1 (u32_16): Invalid bit stride (i)");
    static_assert(BSu32_16::c_D0inWszb == 16,
                  "Set1 (u32_16): Invalid D0inW bit size");
    static_assert(BSu32_16::c_stepD0inWb == stc::Vector<unsigned, 3>{16, 1, 1},
                  "Set1 (u32_16): Invalid D0inW bit steps (u)");
    static_assert(BSu32_16::c_iix_stepD0inWb == stc::Vector<int, 3>{16, 1, 1},
                  "Set1 (u32_16): Invalid D0inW bit steps (i)");
    static_assert(BSu32_16::c_maskD0inW == 0xFFFF,
                  "Set1 (u32_16): Invalid D0inW mask");
    // unsigned (size 8 D0)
    using BSu32_8 = MDBitSet<unsigned, int, unsigned, 8, 8, 16, 16>;
    static_assert(BSu32_8::rank == 3, "Set1 (u32_8): Invalid rank");
    static_assert(BSu32_8::Nb == 2048, "Set1 (u32_8): Invalid num bits");
    static_assert(BSu32_8::WszB == 4,
                  "Set1 (u32_8): Invalid word size (B)");
    static_assert(BSu32_8::Wszb == 32,
                  "Set1 (u32_8): Invalid word size (b)");
    static_assert(BSu32_8::NW == 64, "Set1 (u32_8): Invalid num words");
    static_assert(BSu32_8::Nbrem == 0,
                  "Set1 (u32_8): Invalid remainder bits");
    static_assert(BSu32_8::Nbsrp == 0,
                  "Set1 (u32_8): Invalid surplus bits");
    static_assert(BSu32_8::c_dimsb == stc::Vector<unsigned, 3>{8, 16, 16},
                  "Set1 (u32_8): Invalid dimensions (u)");
    static_assert(BSu32_8::c_iix_dimsb == stc::Vector<int, 3>{8, 16, 16},
                  "Set1 (u32_8): Invalid dimensions (i)");
    static_assert(BSu32_8::c_strideb == stc::Vector<unsigned, 3>{1, 8, 8*16},
                  "Set1 (u32_8): Invalid bit stride (u)");
    static_assert(BSu32_8::c_iix_strideb == stc::Vector<int, 3>{1, 8, 8*16},
                  "Set1 (u32_8): Invalid bit stride (i)");
    static_assert(BSu32_8::c_D0inWszb == 8,
                  "Set1 (u32_8): Invalid D0inW bit size");
    static_assert(BSu32_8::c_stepD0inWb == stc::Vector<unsigned, 3>{8, 1, 1},
                  "Set1 (u32_8): Invalid D0inW bit steps (u)");
    static_assert(BSu32_8::c_iix_stepD0inWb == stc::Vector<int, 3>{8, 1, 1},
                  "Set1 (u32_8): Invalid D0inW bit steps (i)");
    static_assert(BSu32_8::c_maskD0inW == 0xFF,
                  "Set1 (u32_8): Invalid D0inW mask");
    // unsigned (size 24 D0 and size 8 D0inW)
    using BSu32_24a = MDBitSet<unsigned, int, unsigned, 8, 24, 16, 16>;
    static_assert(BSu32_24a::rank == 3, "Set1 (u32_24a): Invalid rank");
    static_assert(BSu32_24a::Nb == 6144, "Set1 (u32_24a): Invalid num bits");
    static_assert(BSu32_24a::WszB == 4,
                  "Set1 (u32_24a): Invalid word size (B)");
    static_assert(BSu32_24a::Wszb == 32,
                  "Set1 (u32_24a): Invalid word size (b)");
    static_assert(BSu32_24a::NW == 192, "Set1 (u32_24a): Invalid num words");
    static_assert(BSu32_24a::Nbrem == 0,
                  "Set1 (u32_24a): Invalid remainder bits");
    static_assert(BSu32_24a::Nbsrp == 0,
                  "Set1 (u32_24a): Invalid surplus bits");
    static_assert(BSu32_24a::c_dimsb == stc::Vector<unsigned, 3>{24, 16, 16},
                  "Set1 (u32_24a): Invalid dimensions (u)");
    static_assert(BSu32_24a::c_iix_dimsb == stc::Vector<int, 3>{24, 16, 16},
                  "Set1 (u32_24a): Invalid dimensions (i)");
    static_assert(BSu32_24a::c_strideb == stc::Vector<unsigned,3>{1, 24, 24*16},
                  "Set1 (u32_24a): Invalid bit stride (u)");
    static_assert(BSu32_24a::c_iix_strideb == stc::Vector<int, 3>{1, 24, 24*16},
                  "Set1 (u32_24a): Invalid bit stride (i)");
    static_assert(BSu32_24a::c_D0inWszb == 8,
                  "Set1 (u32_24a): Invalid D0inW bit size");
    static_assert(BSu32_24a::c_stepD0inWb == stc::Vector<unsigned, 3>{8, 1, 1},
                  "Set1 (u32_24a): Invalid D0inW bit steps (u)");
    static_assert(BSu32_24a::c_iix_stepD0inWb == stc::Vector<int, 3>{8, 1, 1},
                  "Set1 (u32_24a): Invalid D0inW bit steps (i)");
    static_assert(BSu32_24a::c_maskD0inW == 0xFF,
                  "Set1 (u32_24a): Invalid D0inW mask");
    // unsigned (size 24 D0 and size 8 D0inW with remainder)
    using BSu32_24b = MDBitSet<unsigned, int, unsigned, 8, 24, 15, 15>;
    static_assert(BSu32_24b::rank == 3, "Set1 (u32_24b): Invalid rank");
    static_assert(BSu32_24b::Nb == 5400, "Set1 (u32_24b): Invalid num bits");
    static_assert(BSu32_24b::WszB == 4,
                  "Set1 (u32_24b): Invalid word size (B)");
    static_assert(BSu32_24b::Wszb == 32,
                  "Set1 (u32_24b): Invalid word size (b)");
    static_assert(BSu32_24b::NW == 169, "Set1 (u32_24b): Invalid num words");
    static_assert(BSu32_24b::Nbrem == 24,
                  "Set1 (u32_24b): Invalid remainder bits");
    static_assert(BSu32_24b::Nbsrp == 8,
                  "Set1 (u32_24b): Invalid surplus bits");
    static_assert(BSu32_24b::c_dimsb == stc::Vector<unsigned, 3>{24, 15, 15},
                  "Set1 (u32_24b): Invalid dimensions (u)");
    static_assert(BSu32_24b::c_iix_dimsb == stc::Vector<int, 3>{24, 15, 15},
                  "Set1 (u32_24b): Invalid dimensions (i)");
    static_assert(BSu32_24b::c_strideb == stc::Vector<unsigned,3>{1, 24, 24*15},
                  "Set1 (u32_24b): Invalid bit stride (u)");
    static_assert(BSu32_24b::c_iix_strideb == stc::Vector<int, 3>{1, 24, 24*15},
                  "Set1 (u32_24b): Invalid bit stride (i)");
    static_assert(BSu32_24b::c_D0inWszb == 8,
                  "Set1 (u32_24b): Invalid D0inW bit size");
    static_assert(BSu32_24b::c_stepD0inWb == stc::Vector<unsigned, 3>{8, 1, 1},
                  "Set1 (u32_24b): Invalid D0inW bit steps (u)");
    static_assert(BSu32_24b::c_iix_stepD0inWb == stc::Vector<int, 3>{8, 1, 1},
                  "Set1 (u32_24b): Invalid D0inW bit steps (i)");
    static_assert(BSu32_24b::c_maskD0inW == 0xFF,
                  "Set1 (u32_24b): Invalid D0inW mask");
    // unsigned long long (size 16 D0)
    using BSu64_16 = MDBitSet<unsigned long long,
                               long long,
                               unsigned long long,
                               16, 16, 16, 16>;
    static_assert(BSu64_16::rank == 3, "Set1 (u64_16): Invalid rank");
    static_assert(BSu64_16::Nb == 4096, "Set1 (u64_16): Invalid num bits");
    static_assert(BSu64_16::WszB == 8,
                  "Set1 (u64_16): Invalid word size (B)");
    static_assert(BSu64_16::Wszb == 64,
                  "Set1 (u64_16): Invalid word size (b)");
    static_assert(BSu64_16::NW == 64, "Set1 (u64_16): Invalid num words");
    static_assert(BSu64_16::Nbrem == 0,
                  "Set1 (u64_16): Invalid remainder bits");
    static_assert(BSu64_16::Nbsrp == 0,
                  "Set1 (u64_16): Invalid surplus bits");
    static_assert(BSu64_16::c_dimsb ==
                  stc::Vector<unsigned long long, 3>{16, 16, 16},
                  "Set1 (u64_16): Invalid dimensions (u)");
    static_assert(BSu64_16::c_iix_dimsb ==
                  stc::Vector<long long, 3>{16, 16, 16},
                  "Set1 (u64_16): Invalid dimensions (i)");
    static_assert(BSu64_16::c_strideb ==
                  stc::Vector<unsigned long long, 3>{1, 16, 16*16},
                  "Set1 (u64_16): Invalid bit stride (u)");
    static_assert(BSu64_16::c_iix_strideb ==
                  stc::Vector<long long, 3>{1, 16, 16*16},
                  "Set1 (u64_16): Invalid bit stride (i)");
    static_assert(BSu64_16::c_D0inWszb == 16,
                  "Set1 (u64_16): Invalid D0inW bit size");
    static_assert(BSu64_16::c_stepD0inWb ==
                  stc::Vector<unsigned long long, 3>{16, 1, 1},
                  "Set1 (u64_16): Invalid D0inW bit steps (u)");
    static_assert(BSu64_16::c_iix_stepD0inWb ==
                  stc::Vector<long long, 3>{16, 1, 1},
                  "Set1 (u64_16): Invalid D0inW bit steps (i)");
    static_assert(BSu64_16::c_maskD0inW == 0xFFFFull,
                  "Set1 (u64_16): Invalid D0inW mask");
    // unsigned long long (size 24 D0 and size 8 D0inW with remainder)
    using BSu64_24b = MDBitSet<unsigned long long,
                              long long,
                              unsigned long long,
                              8, 24, 15, 15>;
    static_assert(BSu64_24b::rank == 3, "Set1 (u64_24b): Invalid rank");
    static_assert(BSu64_24b::Nb == 5400, "Set1 (u64_24b): Invalid num bits");
    static_assert(BSu64_24b::WszB == 8,
                  "Set1 (u64_24b): Invalid word size (B)");
    static_assert(BSu64_24b::Wszb == 64,
                  "Set1 (u64_24b): Invalid word size (b)");
    static_assert(BSu64_24b::NW == 85, "Set1 (u64_24b): Invalid num words");
    static_assert(BSu64_24b::Nbrem == 24,
                  "Set1 (u64_24b): Invalid remainder bits");
    static_assert(BSu64_24b::Nbsrp == 40,
                  "Set1 (u64_24b): Invalid surplus bits");
    static_assert(BSu64_24b::c_dimsb ==
                  stc::Vector<unsigned long long, 3>{24, 15, 15},
                  "Set1 (u64_24b): Invalid dimensions (u)");
    static_assert(BSu64_24b::c_iix_dimsb ==
                  stc::Vector<long long, 3>{24, 15, 15},
                  "Set1 (u64_24b): Invalid dimensions (i)");
    static_assert(BSu64_24b::c_strideb ==
                  stc::Vector<unsigned long long, 3>{1, 24, 24*15},
                  "Set1 (u64_24b): Invalid bit stride (u)");
    static_assert(BSu64_24b::c_iix_strideb ==
                  stc::Vector<long long, 3>{1, 24, 24*15},
                  "Set1 (u64_24b): Invalid bit stride (i)");
    static_assert(BSu64_24b::c_D0inWszb == 8,
                  "Set1 (u64_24b): Invalid D0inW bit size");
    static_assert(BSu64_24b::c_stepD0inWb ==
                  stc::Vector<unsigned long long, 3>{8, 1, 1},
                  "Set1 (u64_24b): Invalid D0inW bit steps (u)");
    static_assert(BSu64_24b::c_iix_stepD0inWb ==
                  stc::Vector<long long, 3>{8, 1, 1},
                  "Set1 (u64_24b): Invalid D0inW bit steps (i)");
    static_assert(BSu64_24b::c_maskD0inW == 0xFFull,
                  "Set1 (u64_24b): Invalid D0inW mask");
    // unsigned char (size 16 D0).  D0inW size is overruled.
    using BSu8_16 = MDBitSet<unsigned char,
                             int,
                             unsigned,
                             16, 16, 8, 8>;
    static_assert(BSu8_16::rank == 3, "Set1 (u8_16): Invalid rank");
    static_assert(BSu8_16::Nb == 1024, "Set1 (u8_16): Invalid num bits");
    static_assert(BSu8_16::WszB == 1,
                  "Set1 (u8_16): Invalid word size (B)");
    static_assert(BSu8_16::Wszb == 8,
                  "Set1 (u8_16): Invalid word size (b)");
    static_assert(BSu8_16::NW == 128, "Set1 (u8_16): Invalid num words");
    static_assert(BSu8_16::Nbrem == 0,
                  "Set1 (u8_16): Invalid remainder bits");
    static_assert(BSu8_16::Nbsrp == 0,
                  "Set1 (u8_16): Invalid surplus bits");
    static_assert(BSu8_16::c_dimsb ==
                  stc::Vector<unsigned, 3>{16, 8, 8},
                  "Set1 (u8_16): Invalid dimensions (u)");
    static_assert(BSu8_16::c_iix_dimsb ==
                  stc::Vector<int, 3>{16, 8, 8},
                  "Set1 (u8_16): Invalid dimensions (i)");
    static_assert(BSu8_16::c_strideb ==
                  stc::Vector<unsigned, 3>{1, 16, 16*8},
                  "Set1 (u8_16): Invalid bit stride (u)");
    static_assert(BSu8_16::c_iix_strideb ==
                  stc::Vector<int, 3>{1, 16, 16*8},
                  "Set1 (u8_16): Invalid bit stride (i)");
    static_assert(BSu8_16::c_D0inWszb == 8,
                  "Set1 (u8_16): Invalid D0inW bit size");
    static_assert(BSu8_16::c_stepD0inWb ==
                  stc::Vector<unsigned, 3>{8, 1, 1},
                  "Set1 (u8_16): Invalid D0inW bit steps (u)");
    static_assert(BSu8_16::c_iix_stepD0inWb ==
                  stc::Vector<int, 3>{8, 1, 1},
                  "Set1 (u8_16): Invalid D0inW bit steps (i)");
    static_assert(BSu8_16::c_maskD0inW == 0xFF,
                  "Set1 (u8_16): Invalid D0inW mask");
  }

//--Set 2 (set 1st and last bit in set, last in word and last in D0inW)
//  ---W0--- ---W1--- ... --WNW-1-
//  -D0-
//  10011001 10000000 ... 00000001

  {
    int stat2 = 0;
    {
      using BSu32_16 = MDBitSet<unsigned, int, unsigned, 16, 16, 16, 16>;
      BSu32_16 bs;
      bs.setb(BSu32_16::IIxVec{ 0,  0,  0});
      bs.setb(BSu32_16::IIxVec{15,  0,  0});
      bs.setb(BSu32_16::IIxVec{ 0,  1,  0});
      bs.setb(BSu32_16::IIxVec{15,  1,  0});
      bs.setb(BSu32_16::IIxVec{ 0,  2,  0});
      bs.setb(BSu32_16::IIxVec{15, 15, 15});
      if (bs.count() != 6) ++stat2;
      // First W
      if (bs.getW(0) != (((unsigned)1 << 31) |
                         ((unsigned)1 << 16) |
                         ((unsigned)1 << 15) |
                         ((unsigned)1 <<  0))) ++stat2;
      // Second W
      if (bs.getW(1) != ((unsigned)1 <<  0)) ++stat2;
      // Interior W
      for (int cW = 2; cW != 127; ++cW)
        {
          if (bs.getW(cW) != 0) ++stat2;
        }
      // Last W
      if (bs.getW(127) != ((unsigned)1 <<  31)) ++stat2;
      // First D0inW
      if (bs.loadD0inW( 0) != (((unsigned)1 << 15) |
                               ((unsigned)1 <<  0))) ++stat2;
      if (bs.loadD0inW(15) != (((unsigned)1 << 15) |
                               ((unsigned)1 <<  0))) ++stat2;
      // Second D0inW (still first word)
      if (bs.loadD0inW(16) != (((unsigned)1 << 15) |
                               ((unsigned)1 <<  0))) ++stat2;
      if (bs.loadD0inW(31) != (((unsigned)1 << 15) |
                               ((unsigned)1 <<  0))) ++stat2;
      // Third D0inW (second word)
      if (bs.loadD0inW(32) != ((unsigned)1 <<  0)) ++stat2;
      if (bs.loadD0inW(47) != ((unsigned)1 <<  0)) ++stat2;
      // Second last D0inW (last word)
      if (bs.loadD0inW(4064) != (unsigned)0) ++stat2;
      if (bs.loadD0inW(4079) != (unsigned)0) ++stat2;
      // Last D0inW (last word)
      if (bs.loadD0inW(4080) != ((unsigned)1 << 15)) ++stat2;
      if (bs.loadD0inW(4095) != ((unsigned)1 << 15)) ++stat2;
      // Iteration
      int checks = 6;
      int cnt = 0;
      bs.forEachTrueBit(BSu32_16::IIxVec(0),
                        [&checks, &cnt]
                        (BSu32_16::IIxVec a_iv)
                          {
                            ++cnt;
                            if (a_iv == BSu32_16::IIxVec{ 0,  0,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{15,  0,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{ 0,  1,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{15,  1,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{ 0,  2,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{15, 15, 15}) --checks;
                          });
      stat2 += checks;
      if (cnt != 6) ++stat2;
      // Range restricted
      checks = 2;
      cnt = 0;
      bs.forEachTrueBit(BSu32_16::IIxVec(0),
                        BSu32_16::IIxVec{0, 1, 0}, BSu32_16::IIxVec{7, 15, 15},
                        [&checks, &cnt]
                        (BSu32_16::IIxVec a_iv)
                          {
                            ++cnt;
                            // if (a_iv == BSu32_16::IIxVec{ 0,  0,  0}) --checks;
                            // if (a_iv == BSu32_16::IIxVec{15,  0,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{ 0,  1,  0}) --checks;
                            // if (a_iv == BSu32_16::IIxVec{15,  1,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{ 0,  2,  0}) --checks;
                            // if (a_iv == BSu32_16::IIxVec{15, 15, 15}) --checks;
                          });
      stat2 += checks;
      if (cnt != 2) ++stat2;
    }
    if (stat2 != 0)
      {
        pout() << "Failure in set 2" << endl;
        status += stat2;
      }
  }

//--Set 3 (test loading and storing of words)

  {
    int stat3 = 0;
    {
      using BSu32_32 = MDBitSet<unsigned, int, unsigned, 32, 32, 16, 16>;
      BSu32_32 bs;
      bs.storeD0inW(((unsigned)1 << 15) | ((unsigned)1 <<  0), 32);
      int checks = 2;
      int cnt = 0;
      bs.forEachTrueBit(BSu32_32::IIxVec(0),
                        [&checks, &cnt]
                        (BSu32_32::IIxVec a_iv)
                          {
                            ++cnt;
                            if (a_iv == BSu32_32::IIxVec{ 0,  1,  0}) --checks;
                            if (a_iv == BSu32_32::IIxVec{15,  1,  0}) --checks;
                          });
      stat3 += checks;
      if (cnt != 2) ++stat3;
      bs.orD0inW(((unsigned)1 << 14) | ((unsigned)1 <<  1), 32);
      checks = 4;
      cnt = 0;
      bs.forEachTrueBit(BSu32_32::IIxVec(0),
                        [&checks, &cnt]
                        (BSu32_32::IIxVec a_iv)
                          {
                            ++cnt;
                            if (a_iv == BSu32_32::IIxVec{ 0,  1,  0}) --checks;
                            if (a_iv == BSu32_32::IIxVec{ 1,  1,  0}) --checks;
                            if (a_iv == BSu32_32::IIxVec{14,  1,  0}) --checks;
                            if (a_iv == BSu32_32::IIxVec{15,  1,  0}) --checks;
                          });
      stat3 += checks;
      if (cnt != 4) ++stat3;
      bs.storeD0inW(((unsigned)1 << 14) | ((unsigned)1 <<  1), 32);
      checks = 2;
      cnt = 0;
      bs.forEachTrueBit(BSu32_32::IIxVec(0),
                        [&checks, &cnt]
                        (BSu32_32::IIxVec a_iv)
                          {
                            ++cnt;
                            if (a_iv == BSu32_32::IIxVec{ 1,  1,  0}) --checks;
                            if (a_iv == BSu32_32::IIxVec{14,  1,  0}) --checks;
                          });
      stat3 += checks;
      if (cnt != 2) ++stat3;
    }
    {
      using BSu32_16 = MDBitSet<unsigned, int, unsigned, 16, 16, 16, 16>;
      BSu32_16 bs;
      bs.storeD0inW(((unsigned)1 << 15) | ((unsigned)1 <<  0), 16);
      int checks = 2;
      int cnt = 0;
      bs.forEachTrueBit(BSu32_16::IIxVec(0),
                        [&checks, &cnt]
                        (BSu32_16::IIxVec a_iv)
                          {
                            ++cnt;
                            if (a_iv == BSu32_16::IIxVec{ 0,  1,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{15,  1,  0}) --checks;
                          });
      stat3 += checks;
      if (cnt != 2) ++stat3;
      bs.orD0inW(((unsigned)1 << 14) | ((unsigned)1 <<  1), 16);
      checks = 4;
      cnt = 0;
      bs.forEachTrueBit(BSu32_16::IIxVec(0),
                        [&checks, &cnt]
                        (BSu32_16::IIxVec a_iv)
                          {
                            ++cnt;
                            if (a_iv == BSu32_16::IIxVec{ 0,  1,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{ 1,  1,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{14,  1,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{15,  1,  0}) --checks;
                          });
      stat3 += checks;
      if (cnt != 4) ++stat3;
      bs.storeD0inW(((unsigned)1 << 14) | ((unsigned)1 <<  1), 16);
      checks = 2;
      cnt = 0;
      bs.forEachTrueBit(BSu32_16::IIxVec(0),
                        [&checks, &cnt]
                        (BSu32_16::IIxVec a_iv)
                          {
                            ++cnt;
                            if (a_iv == BSu32_16::IIxVec{ 1,  1,  0}) --checks;
                            if (a_iv == BSu32_16::IIxVec{14,  1,  0}) --checks;
                          });
      stat3 += checks;
      if (cnt != 2) ++stat3;
    }
    if (stat3 != 0)
      {
        pout() << "Failure in set 3" << endl;
        status += stat3;
      }
  }

//--Set 4 (test accessing words in a loop)

  {
    int stat4 = 0;
    {
      using BSu32_16 = MDBitSet<unsigned, int, unsigned, 16, 16, 16, 16>;
      BSu32_16 bs;
      bs.setAll();
      int checks = bs.c_dimsWD0inW.product();
      bs.forEachCWord(BSu32_16::IIxVec(0),
                     [&checks]
                     (const BSu32_16::IIxVec& a_iv,
                      BSu32_16::word_type     a_mask,
                      BSu32_16::word_type     a_word)
                       {
                         if (a_mask == a_word) --checks;
                       });
      stat4 += checks;
      checks = bs.c_dimsWD0inW.product();
      bs.forEachWord(BSu32_16::IIxVec(0),
                     [&checks]
                     (const BSu32_16::IIxVec& a_iv,
                      BSu32_16::word_type     a_mask,
                      BSu32_16::word_type&    a_word)
                       {
                         if (a_iv[0] == 0) --checks;
                         //       Set 1 bit   Should be ignored
                         a_word = (1<<1)    | (1<<16);
                       });
      stat4 += checks;
      checks = bs.c_dimsWD0inW.product();
      int cnt = 0;
      bs.forEachTrueBit(BSu32_16::IIxVec(0),
                        [&checks, &cnt]
                        (BSu32_16::IIxVec a_iv)
                          {
                            ++cnt;
                            if (a_iv[0] == 1) --checks;
                          });
      stat4 += checks;
      if (cnt != bs.c_dimsWD0inW.product()) ++stat4;
      if (bs.count() != bs.c_dimsWD0inW.product()) ++stat4;
      // With range
      bs.setAll();
      checks = 16;
      bs.forEachWord(BSu32_16::IIxVec(0),
                     BSu32_16::IIxVec{ 9, 1, 12},
                     BSu32_16::IIxVec{14, 4, 15},
                     [&checks]
                     (const BSu32_16::IIxVec& a_iv,
                      BSu32_16::word_type     a_mask,
                      BSu32_16::word_type&    a_word)
                       {
                         if (a_mask == (BSu32_16::word_type)0x7E00 &&
                             a_word == (BSu32_16::word_type)0xFFFF &&
                             a_iv[0] == 0 &&
                             a_iv[1] >=  1 && a_iv[1] <= 4 &&
                             a_iv[2] >= 12 && a_iv[1] <= 15) --checks;
                         a_word = (BSu32_16::word_type)0;
                       });
      stat4 += checks;
      // You can only modify the bits in range even though you are given all
      // bits.  So where we set a_word = 0 in the above, only the bits 9--14
      // should be zeroed.
      checks = 16;
      bs.forEachWord(BSu32_16::IIxVec(0),
                     BSu32_16::IIxVec{ 9, 1, 12},
                     BSu32_16::IIxVec{14, 4, 15},
                     [&checks]
                     (const BSu32_16::IIxVec& a_iv,
                      BSu32_16::word_type     a_mask,
                      BSu32_16::word_type     a_word)
                       {
                         if (a_word == (BSu32_16::word_type)0x81FF) --checks;
                       });
      stat4 += checks;
      checks = bs.c_dimsWD0inW.product();
      bs.forEachWord(BSu32_16::IIxVec(0),
                     [&checks]
                     (const BSu32_16::IIxVec& a_iv,
                      BSu32_16::word_type     a_mask,
                      BSu32_16::word_type     a_word)
                       {
                         if (a_iv >= BSu32_16::IIxVec{0, 1, 12} &&
                             a_iv <= BSu32_16::IIxVec{0, 4, 15})
                           {
                             if (a_word == (BSu32_16::word_type)0x81FF)
                               --checks;
                           }
                         else
                           {
                             if (a_word == (BSu32_16::word_type)0xFFFF)
                               --checks;
                           }
                       });
      stat4 += checks;

    }
    if (stat4 != 0)
      {
        pout() << "Failure in set 4" << endl;
        status += stat4;
      }
  }

//--Set 5 (minbox)

  {
    int stat5 = 0;
    {
      using BSu32_16 = MDBitSet<unsigned, int, unsigned, 16, 16, 16, 16>;
      BSu32_16 bs;
      BSu32_16::IIxVec lb, ub;
      // No bits
      bs.minBox(lb, ub);
      if (lb != BSu32_16::c_iix_dimsb) ++stat5;
      if (ub != BSu32_16::IIxVec(-1)) ++stat5;
      // {0, 0, 0}
      bs.setb(BSu32_16::IIxVec(0));
      bs.minBox(lb, ub);
      if (lb != BSu32_16::IIxVec(0)) ++stat5;
      if (ub != BSu32_16::IIxVec(0)) ++stat5;
      // All bits
      bs.setAll();
      bs.minBox(lb, ub);
      if (lb != BSu32_16::IIxVec(0)) ++stat5;
      if (ub != BSu32_16::c_iix_dimsb - 1) ++stat5;
      // Various...
      bs.clearAll();
      bs.setb(BSu32_16::IIxVec{0, 1, 0});
      bs.minBox(lb, ub);
      if (lb != BSu32_16::IIxVec{0, 1, 0}) ++stat5;
      if (ub != BSu32_16::IIxVec{0, 1, 0}) ++stat5;
      bs.setb(BSu32_16::IIxVec{0, 2, 15});
      bs.minBox(lb, ub);
      if (lb != BSu32_16::IIxVec{0, 1, 0}) ++stat5;
      if (ub != BSu32_16::IIxVec{0, 2, 15}) ++stat5;
      bs.setb(BSu32_16::IIxVec{1, 1, 0});
      bs.minBox(lb, ub);
      if (lb != BSu32_16::IIxVec{0, 1, 0}) ++stat5;
      if (ub != BSu32_16::IIxVec{1, 2, 15}) ++stat5;
      bs.clearAll();
      bs.setb(BSu32_16::IIxVec{0, 0, 7});
      bs.setb(BSu32_16::IIxVec{0, 0, 8});
      bs.setb(BSu32_16::IIxVec{1, 1, 7});
      bs.setb(BSu32_16::IIxVec{1, 1, 8});
      bs.minBox(lb, ub);
      if (lb != BSu32_16::IIxVec{0, 0, 7}) ++stat5;
      if (ub != BSu32_16::IIxVec{1, 1, 8}) ++stat5;
    }
    if (stat5 != 0)
      {
        pout() << "Failure in set 5" << endl;
        status += stat5;
      }
  }
  
//--Set 5 (three-way)

  {
    int stat5 = 0;
    {
      using BSu32_16 = MDBitSet<unsigned, int, unsigned, 16, 16, 16, 16>;
      BSu32_16 bs;
      auto bs3w = bs.makeThreeWay();
      std::cout << bs3w.c_dimsb << std::endl;
      std::cout << bs3w.c_D0inWszb << std::endl;
    }
    if (stat5 != 0)
      {
        pout() << "Failure in set 5" << endl;
        status += stat5;
      }
  }

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
