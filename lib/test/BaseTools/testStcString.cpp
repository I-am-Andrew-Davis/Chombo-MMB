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
#include <sstream>

#if CXXSTD>=14
#define D_DECL6(a,b,c,d,e,f) a,b,c
#include "StcString.H"
#include "CH_Hash.H"
#endif
#include "parstream.H"
#include "SPMD.H"

#include "UsingBaseNamespace.H"

/// Prototypes:

void
parseTestOptions(int argc, char* argv[]);

int
testStcString();

/// Global variables for handling output:
static const char *pgmname = "testStcString" ;
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
  int ret = testStcString();

  if (ret == 0)
    {
      if (verbose)
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

/*------------------------------------------------------------------------------
 * Tests
 *----------------------------------------------------------------------------*/

int testStcString()
{
  int status = 0;
#if CXXSTD>=14

//--Tests

  // Set 1 (example)
  if (verbose)
    {
      stc::String<8> strA("abcd");
      strA += "efghik";
      std::cout << "A: " << strA << std::endl;
    }

  // Set 2
  // Test construction and equivalence using c-style strings
  {
    int lstat = 0;
    stc::String<4> strA("abcd");
    if (strA.size() != 5) ++lstat;
    if (strA.length() != 4) ++lstat;
    if (strA[0] != 'a') ++lstat;
    if (strA[1] != 'b') ++lstat;
    if (strA[2] != 'c') ++lstat;
    if (strA[3] != 'd') ++lstat;
    if (strA[4] != '\0') ++lstat;
    if (!(strA == "abcd")) ++lstat;
    if (!(strA == "abcde")) ++lstat;
    if   (strA == "abc") ++lstat;
    if   (strA != "abcd") ++lstat;
    if   (strA != "abcde") ++lstat;
    if (!(strA != "abc")) ++lstat;
    if (!("abcd"  == strA)) ++lstat;
    if (!("abcde" == strA)) ++lstat;
    if   ("abc"   == strA) ++lstat;
    if   ("abcd"  != strA) ++lstat;
    if   ("abcde" != strA) ++lstat;
    if (!("abc"   != strA)) ++lstat;
    strA = "xyz";
    if (strA.size() != 5) ++lstat;
    if (strA.length() != 3) ++lstat;
    if (strA[0] != 'x') ++lstat;
    if (strA[1] != 'y') ++lstat;
    if (strA[2] != 'z') ++lstat;
    if (strA[3] != '\0') ++lstat;
    if (strA[4] != '\0') ++lstat;
    if (!(strA == "xyz")) ++lstat;
    if (  strA == "xyza") ++lstat;
    if   (strA == "xy") ++lstat;
    if   (strA != "xyz") ++lstat;
    if (!(strA != "xyza")) ++lstat;
    if (!(strA != "xy")) ++lstat;
    if (!("xyz"  == strA)) ++lstat;
    if   ("xyza" == strA) ++lstat;
    if   ("xy"   == strA) ++lstat;
    if   ("xyz"  != strA) ++lstat;
    if (!("xyza" != strA)) ++lstat;
    if (!("xy"   != strA)) ++lstat;
    // Test a null str
    strA = "";
    if (strA.length() != 0) ++lstat;
    if (strA[0] != '\0') ++lstat;
    if (strA != "") ++lstat;
    stc::String<4> strB("");
    if (strB.length() != 0) ++lstat;
    if (strB[0] != '\0') ++lstat;
    if (strB != "") ++lstat;
    if (lstat != 0)
      {
        pout() << "Failure in set 2" << std::endl;
        status += lstat;
      }
  }

  // Set 3
  // Test construction and equivalence using std::string
  {
    int lstat = 0;
    stc::String<4> strA(std::string("abcd"));
    if (strA.size() != 5) ++lstat;
    if (strA.length() != 4) ++lstat;
    if (!(strA == "abcd")) ++lstat;
    if (!(strA == std::string("abcd"))) ++lstat;
    if (!(strA == std::string("abcde"))) ++lstat;
    if   (strA == std::string("abc")) ++lstat;
    if   (strA != std::string("abcd")) ++lstat;
    if   (strA != std::string("abcde")) ++lstat;
    if (!(strA != std::string("abc"))) ++lstat;
    if (!(std::string("abcd")  == strA)) ++lstat;
    if (!(std::string("abcde") == strA)) ++lstat;
    if   (std::string("abc")   == strA) ++lstat;
    if   (std::string("abcd")  != strA) ++lstat;
    if   (std::string("abcde") != strA) ++lstat;
    if (!(std::string("abc")   != strA)) ++lstat;
    if (lstat != 0)
      {
        pout() << "Failure in set 3" << std::endl;
        status += lstat;
      }
  }

  // Set 4
  // Constexpr
  {
    constexpr stc::String<4> strA('b');
    static_assert(strA[0] == 'b', "Failed constexpr construction (4)");
    constexpr stc::String<4> strB("abc");
    static_assert(strB == "abc", "Failed constexpr construction (4)");
    constexpr stc::String<4> strC(strB);
    static_assert(strC == "abc", "Failed constexpr construction (4)");
  }

  // Set 5
  // += concatentation
  {
    int lstat = 0;
    stc::String<12> strA("abc");
    constexpr stc::String<12> strB("ghi");
    constexpr stc::String<12> strC("abcdefghijkl");
    strA += 'd';
    if (strA != "abcd") ++lstat;
    strA += "ef";
    if (strA != "abcdef") ++lstat;
    strA += strB;
    if (strA != "abcdefghi") ++lstat;
    strA += std::string("jkl");
    if (strA != "abcdefghijkl") ++lstat;
    strA = "abc";
    strA += "defghijklmnopqrstuvwxy";
    if (strA != strC) ++lstat;
    strA = "abc";
    strA += std::string("defghijklmnopqrstuvwxy");
    if (strA != strC) ++lstat;
    strA += "xyz";
    if (strA != strC) ++lstat;
    strA += std::string("xyz");
    if (strA != strC) ++lstat;
    strA += 'z';
    if (strA != strC) ++lstat;
    strA += strB;
    if (strA != strC) ++lstat;
    if (lstat != 0)
      {
        pout() << "Failure in set 5" << std::endl;
        status += lstat;
      }
  }

  // Set 6
  // + concatentation
  {
    int lstat = 0;
    constexpr stc::String<4> strA("abc");
    constexpr stc::String<8> strB("defghij");
    const auto strC = strA + strB;  // Generates C with size max(4, 8)
    if (strC != "abcdefgh") ++lstat;
    if (strC.size() != 9) ++lstat;
    if (strC.length() != 8) ++lstat;
    // While strD can store 8 chars, the temporary from operator + only has
    // storage for 4 chars
    stc::String<8> strD = strA + "defj";
    if (strD.length() != 4) ++lstat;
    if (strD == stc::String<8>("abcdefjh")) ++lstat;
    if (strD != stc::String<8>("abcd")) ++lstat;
    strD = strA + 'd';
    if (strD.length() != 4) ++lstat;
    if (strD != stc::String<8>("abcd")) ++lstat;
    strD = 'd' + strA;
    if (strD.length() != 4) ++lstat;
    if (strD != stc::String<8>("dabc")) ++lstat;
    constexpr stc::String<12> strE("abc");
    strD = strE + "def";
    if (strD.length() != 6) ++lstat;
    if (strD != stc::String<8>("abcdef")) ++lstat;
    strD = "defghi" + strE;
    if (strD.length() != 8) ++lstat;
    if (strD != stc::String<8>("defghiab")) ++lstat;
    strD = strE + std::string("def");
    if (strD.length() != 6) ++lstat;
    if (strD != stc::String<8>("abcdef")) ++lstat;
    strD = std::string("defghi") + strE;
    if (strD.length() != 8) ++lstat;
    if (strD != stc::String<8>("defghiab")) ++lstat;
    if (lstat != 0)
      {
        pout() << "Failure in set 6" << std::endl;
        status += lstat;
      }
  }

  // Set 7
  // Constructed aliases
  {
    int lstat = 0;
    constexpr stc::String<8> strA("abcd");
    // An alias is restricted in size up to the extra '\0' character in the
    // source array.  This means if an offset is provided, the size must be
    // reduced.  The alias also must be large enough to store whatever is
    // currently in the source array.
    stc::String<7, stc::VectorConstAliasImpl> strB(strA, 1);
    if (strB != stc::String<3>("bcd")) ++lstat;
    stc::String<8> strC("abc");
    stc::String<7, stc::VectorAliasImpl> strD(strC, 1);
    strC += "def";
    if (strD != "bcdef") ++lstat;
    // std::cout << strK << ' ' << strK.size() << ' ' << strK.length()
    //           << std::endl;
    stc::String<7, stc::VectorConstAliasImpl> strE(strD, 0);
    if (strE != "bcdef") ++lstat;
    strD[0] = 'z';
    if (strC != "azcdef") ++lstat;
    if (strE != "zcdef") ++lstat;
    strD = "morethan7char";
    if (strC != stc::String<8>("amoretha")) ++lstat;
    if (strE != stc::String<7>("moretha")) ++lstat;
    if (lstat != 0)
      {
        pout() << "Failure in set 7" << std::endl;
        status += lstat;
      }
  }

  // Set 8
  // Read input
  {
    int lstat = 0;
    std::istringstream input("tuxedo penguin");
    stc::String<6> strA, strB;
    input >> strA >> strB;
    if (strA != "tuxedo") ++lstat;
    if (strB != stc::String<6>("pengui")) ++lstat;
    if (verbose)
      {
        std::cout << strB << "n wearing " << strA << std::endl;
      }
    if (lstat != 0)
      {
        pout() << "Failure in set 7" << std::endl;
        status += lstat;
      }
  }

  // Set 9
  // Hashes
  {
    int lstat = 0;
    stc::String<8> strA("cons8");
    stc::String<6> strB("cons8");
    stc::String<8, stc::VectorConstAliasImpl> strC(strA);
    auto hashA = CH_Hash::XXHash32<stc::String<8>>{}(strA);
    auto hashB = CH_Hash::hash_XXHash32(strB);
    auto hashC = CH_Hash::hash_XXHash32(strC);
    if (hashA != hashB) ++lstat;
    if (hashA != hashC) ++lstat;
    strA += "toPrim";
    strB += "toPrim";
    hashA = CH_Hash::XXHash32<stc::String<8>>{}(strA);
    hashB = CH_Hash::hash_XXHash32(strB);
    hashC = CH_Hash::hash_XXHash32(strC);
    if (hashA == hashB) ++lstat;
    if (hashA != hashC) ++lstat;
    if (lstat != 0)
      {
        pout() << "Failure in set 8" << std::endl;
        status += lstat;
      }
  }

#endif  /* CXXSTD>=14 */
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
