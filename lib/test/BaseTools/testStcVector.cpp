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

#if CXXSTD>=14
#define D_DECL6(a,b,c,d,e,f) a,b,c
#include "StcVector.H"
#endif
#include "parstream.H"
#include "SPMD.H"

#include "UsingBaseNamespace.H"

/// Prototypes:

void
parseTestOptions(int argc, char* argv[]);

int
testStcVector();

/// Global variables for handling output:
static const char *pgmname = "testStcVector" ;
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
  int ret = testStcVector();

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

#if CXXSTD>=14
// Dummy IntVect class with SpaceDim=3
constexpr int c_SpaceDim = 3;
using IntVect = stc::IVec<c_SpaceDim>;
constexpr auto IntVect_zero = stc::make_IVec<c_SpaceDim>::zero();
constexpr auto IntVect_unit = stc::make_IVec<c_SpaceDim>::unit();
constexpr auto IntVect_basis(const size_t a_dir) noexcept
{
  return stc::make_IVec<c_SpaceDim>::basis(a_dir);
}
template <typename Op>
constexpr auto IntVect_make_array(Op&& a_op) noexcept
{
  return stc::Make_array<int, c_SpaceDim, Op>::with(std::forward<Op>(a_op));
}
#endif

/*------------------------------------------------------------------------------
 * WARNING!  Do not define constexpr values in real code unless you are sure and
 * test well.  These cannot evaluate to constexpr values:
 * Intel (v 19.0.0.117)
 *   IntVect_basis
 *   IntVect_make_array
 *   stc::cross
 *----------------------------------------------------------------------------*/

int testStcVector()
{
  int status = 0;
#if CXXSTD>=14

//--Tests

  // Preliminary assetions
  static_assert(sizeof(IntVect) == c_SpaceDim*sizeof(int),
                "Bad assumption for sizeof(IntVect)");

  // Set (1a)
  {
    constexpr IntVect iv{ 0, 0, -1, 0, 0, 0 };
    if (verbose)
      {
        pout() << "For IntVect: " << iv << std::endl;
      }
    static_assert(iv[2] == -1, "Bad value (1a)");
    constexpr int sum1 = iv.sum();
    constexpr int norm = iv.norm1();
    if (verbose)
      {
        pout() << "  sum   : " << sum1 << std::endl;
        pout() << "  1-norm: " << norm << std::endl;
      }
    static_assert(sum1 == -1, "Bad value (1a)");
    static_assert(norm == 1, "Bad value (1a)");
    constexpr int sum2 = sum(iv);
    static_assert(sum1 == sum2, "Bad value (1a)");
    constexpr IntVect ivB(0, 1, 2);
    static_assert(ivB == IntVect{ 0, 1, 2 }, "Bad value (1a)");
  }
  
  // Set (1b)
  {
    constexpr IntVect iv(3);
    if (verbose)
      {
        pout() << "For IntVect: " << iv << std::endl;
      }
    static_assert(iv == IntVect{ 3, 3, 3, 3, 3, 3 }, "Bad construction (1b)");
  }

  // Set (1c)
  // Extension to larger vector
  {
    constexpr stc::Vector<int, 3> iva{ 1, 2, 3 };
    constexpr stc::Vector<int, 4> ivb(iva, 4);
    static_assert(ivb ==  stc::Vector<int, 4>{ 1, 2, 3, 4 }, "Bad value (1c)");
    constexpr stc::Vector<int, 4> ivc(4, iva);
    static_assert(ivc ==  stc::Vector<int, 4>{ 4, 1, 2, 3 }, "Bad value (1c)");
    constexpr stc::Vector<int, 4> ivd(4, 5, iva);
    static_assert(ivd ==  stc::Vector<int, 4>{ 4, 5, 1, 2 }, "Bad value (1c)");
    constexpr stc::Vector<int, 4> ive(4, 5, 6, iva);
    static_assert(ive ==  stc::Vector<int, 4>{ 4, 5, 6, 1 }, "Bad value (1c)");
  }

  // Set (1d)
  // Extension to larger vector
  {
    constexpr stc::Vector<Real, 3> iva{ 1.5, 2.5, 3.5 };
    constexpr stc::Vector<Real, 4> ivb(iva, 4.5);
    static_assert(ivb ==  stc::Vector<Real, 4>{ 1.5, 2.5, 3.5, 4.5 },
                  "Bad value (1d)");
    constexpr stc::Vector<Real, 4> ivc(4.5, iva);
    static_assert(ivc ==  stc::Vector<Real, 4>{ 4.5, 1.5, 2.5, 3.5 },
                  "Bad value (1d)");
    constexpr stc::Vector<Real, 4> ivd(4.5, 5.5, iva);
    static_assert(ivd ==  stc::Vector<Real, 4>{ 4.5, 5.5, 1.5, 2.5 },
                  "Bad value (1d)");
    constexpr stc::Vector<Real, 4> ive(4.5, 5.5, 6.5, iva);
    static_assert(ive ==  stc::Vector<Real, 4>{ 4.5, 5.5, 6.5, 1.5 },
                  "Bad value (1d)");
  }

  // Set (2)
  {
    constexpr IntVect iva = abs(IntVect{ 0, 0, -1, 0, 0, 0 })*-2;
    static_assert(iva[2] == -2, "Bad value (2)");
    static_assert(stc::abs(IntVect{ 0, 0, -1, 0, 0, 0 })*-2 ==
                  IntVect{ 0, 0, -2, 0, 0 }, "Bad absolute");
    constexpr IntVect ivb = IntVect{0, 0, 1, 0, 0, 0 }*-2;
    static_assert(iva == ivb, "Bad value (2)");
    if (verbose)
      {
        pout() << "For IntVect: " << iva  << " + " << ivb << std::endl;
      }
    constexpr int sum1 = sum(iva + ivb);
    constexpr int norm = norm1(iva + ivb);
    constexpr int dot1 = dot(iva, ivb);
    if (verbose)
      {
        pout() << "  sum   : " << sum1 << std::endl;
        pout() << "  1-norm: " << norm << std::endl;
        pout() << "  dot   : " << dot1 << std::endl;
      }
    static_assert(sum1 == -4, "Bad value (2)");
    static_assert(norm == 4, "Bad value (2)");
    static_assert(dot1 == 4, "Bad value (2)");
  }

  // Set (3)
  {
    constexpr int norm = norm1(IntVect{1, 2, 3, 4});
    if (verbose)
      {
        pout() << "For IntVect: " << IntVect{1, 2, 3, 4} << std::endl;
        pout() << "  1-norm: " << norm << std::endl;
      }
    static_assert(norm == 6, "Bad value (3)");
  }

  // Set (4)
  {
#ifndef __INTEL_COMPILER
    constexpr IntVect iv = 2*IntVect_basis(2);
    static_assert(iv[0] == 0, "Bad value (4)");
    static_assert(iv[1] == 0, "Bad value (4)");
    static_assert(iv[2] == 2, "Bad value (4)");
#else
    IntVect iv = 2*IntVect_basis(2);
    if (iv != IntVect{ 0, 0, 2 }) ++status;
#endif
  }

  // Set (5)
  {
    // This is a representation
    constexpr auto rep = 2*std::move(IntVect_unit);  // This is an expression
    // Expressions are not IntVects but you can still use operator[].  However
    // each application of [] re-evaluates the expression.
    constexpr int i0 = rep[0];
    static_assert(i0 == 2, "Bad value (5)");
    // Note: IntVect_unit is an lvalue.  These propagate as references in the
    // expressions.  Rvalues are copied in the expressions so using move allows
    // constexpr to work.  Yes, you read that right, here move causes copy.
    // This also means that an expression can live past a scope.
    constexpr IntVect iv = rep;                      // Evaluates the expression
    if (verbose)
      {
        pout() << "For IntVect: " << iv << std::endl;
      }
    constexpr IntVect ivb(rep);
    static_assert(iv == ivb, "Bad value (5)");
    constexpr int sum1 = sum(iv);
    constexpr int sum2 = sum(rep);
    static_assert(sum1 == sum2, "Bad value (5)");
    constexpr int prod = product(iv);
    constexpr int norm = norm1(iv);
    if (verbose)
      {
        pout() << "  sum   : " << sum1 << std::endl;
        pout() << "  prod  : " << prod << std::endl;
        pout() << "  1-norm: " << norm << std::endl;
      }
    static_assert(sum1 == 6, "Bad value (5)");
    static_assert(prod == 8, "Bad value (5)");
    static_assert(norm == 6, "Bad value (5)");
  }

  // Set (5a)
  // Modulus
  {
    IntVect iva{ 1, 2, 3 };
    IntVect ivb = iva % 2;
    if (ivb != IntVect{ 1, 0, 1}) ++status;
    ivb = iva % IntVect{ 2, 3, 4 };
    if (ivb != iva) ++status;
    ivb = 2 % iva;
    if (ivb != IntVect{ 0, 0, 2 }) ++status;
  }

  // Set (6)
  // Example: multiply by basis
  {
    IntVect iva{1, 2, 3};
    IntVect ivb(IntVect_zero);
    for (int dir = 0; dir != c_SpaceDim; ++dir)
      {
        ivb += iva*IntVect_basis(dir);  // Note use of rvalue in expression
      }
    if (ivb != iva) ++status;
    if (verbose)
      {
        pout() << "After IV ops: " << ivb << std::endl;
      }
  }

  // Set (7)
  // Example: using IntVect::Unit and IntVect::Zero
  {
    IntVect iva(2*IntVect::Unit - 1);
    IntVect ivb(2*IntVect::Zero);
    if (iva != IntVect_unit) ++status;
    if (ivb != IntVect_zero) ++status;
    if (verbose)
      {
        pout() << "Unit: " << iva << "  " << IntVect::Unit << std::endl;
        pout() << "Zero: " << ivb << "  " << IntVect::Zero << std::endl;
      }
  }

  // Set (8)
  // Min/max testing
  {
    constexpr IntVect iva{ -1, 0, 5, -6 };
    auto exm = stc::min(iva, IntVect_zero);  // Expression
    if (exm != IntVect{ -1, 0, 0, -6 }) ++status;
    IntVect ivx = stc::max(iva, IntVect_zero);
    if (ivx != IntVect{ 0, 0, 5, 0 }) ++status;
    IntVect ivm = stc::min(iva, 0);
    if (ivm != IntVect{ -1, 0, 0, -6 }) ++status;
    ivm = stc::min(0, iva);
    if (ivm != IntVect{ -1, 0, 0, -6 }) ++status;
    ivx = stc::max(iva, 0);
    if (ivx != IntVect{ 0, 0, 5, 0 }) ++status;
    ivx = stc::max(-1, iva);
    if (ivx != iva) ++status;
  }

  // Set (9)
  // Custom makers (for constexpr, need struct)
  {
    struct CustomMake
    {
      constexpr CustomMake(const int a_j)
      :
      m_j(a_j)
        { }
      constexpr int operator()(const size_t a_idx) const noexcept
        {
          return -m_j + a_idx*a_idx;
        }
      const int m_j;
    };
    constexpr IntVect ivcm(IntVect_make_array(CustomMake(2)));
#ifndef __INTEL_COMPILER
    static_assert(ivcm == IntVect{ -2, -1, 2 }, "Bad value");
#endif
    if (verbose)
      {
        pout() << "Custom: " << ivcm << std::endl;
      }
  }

  // Set (10)
  // Comparisons
  {
    constexpr IntVect iva{1, 2, 3};
    // <
    constexpr bool blt1 = IntVect_unit < iva;
    if (blt1) ++status;
    constexpr bool blt2 = IntVect_zero < iva;
    if (!blt2) ++status;
    constexpr bool blt3 = 1 < iva;
    if (blt3) ++status;
    constexpr bool blt4 = 0 < iva;
    if (!blt4) ++status;
    constexpr bool blt5 = iva < 3;
    if (blt5) ++status;
    constexpr bool blt6 = iva < 4;
    if (!blt6) ++status;
    // <=
    constexpr bool ble1 = IntVect_unit <= iva;
    if (!ble1) ++status;
    constexpr bool ble2 = 1 <= iva;
    if (!ble2) ++status;
    constexpr bool ble3 = 2 <= iva;
    if (ble3) ++status;
    constexpr bool ble4 = iva <= 3;
    if (!ble4) ++status;
    constexpr bool ble5 = iva <= 2;
    if (ble5) ++status;
    // >
    static_assert(!(iva > IntVect_unit), "Bad comparison");
    // >=
    static_assert(iva >= IntVect_unit, "Bad comparison");

    static_assert(!iva.lexLT(IntVect_unit), "Bad comparison");
    static_assert(iva.lexGT(IntVect_unit), "Bad comparison");
    static_assert(iva.lexLT(2*IntVect_unit), "Bad comparison");
    static_assert(!iva.lexGT(2*IntVect_unit), "Bad comparison");

    static_assert(lexLT(IntVect{1, 1, 1}, 2*IntVect_unit), "Bad comparison");
    static_assert(lexLT(IntVect{2, 1, 1}, 2*IntVect_unit), "Bad comparison");
    static_assert(lexLT(IntVect{2, 2, 1}, 2*IntVect_unit), "Bad comparison");
    static_assert(!lexLT(IntVect{2, 2, 2}, 2*IntVect_unit), "Bad comparison");
    static_assert(lexGT(IntVect{2, 2, 2}, IntVect_unit), "Bad comparison");
    static_assert(lexGT(IntVect{1, 2, 2}, IntVect_unit), "Bad comparison");
    static_assert(lexGT(IntVect{1, 1, 2}, IntVect_unit), "Bad comparison");
    static_assert(!lexGT(IntVect{1, 1, 1}, IntVect_unit), "Bad comparison");
  }

  // Set (11)
  // Coarsen
  {
    IntVect iv0cr = stc::coarsen(IntVect_zero, 4);
    if (iv0cr != IntVect_zero) ++status;
    constexpr IntVect iva{3, 3, 3};
    if (stc::coarsen(iva, 4) != IntVect_zero) ++status;
    if (stc::coarsen(IntVect{4, 4, 4}, 4) != IntVect_unit) ++status;
    if (stc::coarsen(IntVect{-1, -1, -1}, 4) != -IntVect_unit) ++status;
    if (stc::coarsen(IntVect{-4, -4, -4}, 4) != -IntVect_unit) ++status;
    if (stc::coarsen(IntVect{-5, -5, -5}, 4) != -2*IntVect_unit) ++status;
    if (stc::coarsen(IntVect{-5, -1, 4}, 4) != IntVect{-2, -1, 1}) ++status;
    // Following should fail
    // if (coarsen(IntVect{-5, -1, 4}, 3) != IntVect{-2, -1, 1}) ++status;
    IntVect ivb{-5, -1, 4};
    ivb.coarsen(4);
    if (ivb != IntVect{-2, -1, 1}) ++status;
    // Testing with a vector definition of coarsening factor
    if (stc::coarsen(IntVect{-5, -1, 4}, IntVect{4, 1, 2}) !=
        IntVect{-2, -1, 2}) ++status;
    ivb = {-5, -1, 4};
    if (ivb.coarsen(IntVect{4, 1, 2}) != IntVect{-2, -1, 2}) ++status;

    //--Compile-time
    constexpr IntVect ivc = stc::clog2(IntVect{1, 2, 4});
    static_assert(ivc == IntVect{0, 1, 2}, "Bad log2");
    static_assert(stc::ccoarsen_log2(IntVect{8, 8, 8}, ivc) == IntVect{8, 4, 2},
                  "Bad coarsen");
    static_assert(stc::ccoarsen(IntVect{8, 8, 8}, IntVect{1, 2, 4}) ==
                  IntVect{8, 4, 2}, "Bad coarsen");
    static_assert(stc::ccoarsen(IntVect{-4, 1, 3}, IntVect{4, 2, 1}) ==
                  IntVect{-1, 0, 3}, "Bad coarsen");
    static_assert(stc::ccoarsen(IntVect{-4, 1, 3}, 2) == IntVect{-2, 0, 1},
                  "Bad coarsen");
    static_assert(stc::ccoarsen_log2(IntVect{-4, 1, 3}, 1) == IntVect{-2, 0, 1},
                  "Bad coarsen");
  }

  // Set (12)
  // Cross
  {
#ifndef __INTEL_COMPILER
    static_assert(cross(IntVect_basis(0), IntVect_basis(1)) == IntVect_basis(2),
                  "Bad cross product");
    static_assert(cross(IntVect_basis(0), IntVect_basis(1)) == IntVect_basis(2),
                  "Bad cross product");
    static_assert(cross(IntVect_basis(1), IntVect_basis(2)) == IntVect_basis(0),
                  "Bad cross product");
    static_assert(cross(2*IntVect_basis(2), 2*IntVect_basis(0)) ==
                  4*IntVect_basis(1), "Bad cross product");
    static_assert(cross(IntVect_basis(1), IntVect_basis(0)) ==
                  -IntVect_basis(2), "Bad cross product");
    static_assert(cross(IntVect_basis(2), IntVect_basis(1)) ==
                  -IntVect_basis(0), "Bad cross product");
    static_assert(cross(2*IntVect_basis(0), 2*IntVect_basis(2)) ==
                  -4*IntVect_basis(1), "Bad cross product");
#endif
    constexpr IntVect iva{ 1,  2,  3};
    constexpr IntVect ivb{-1, -2, -3};
#ifndef __INTEL_COMPILER
    constexpr IntVect ivc = cross(iva, ivb);
    static_assert(ivc == IntVect_zero, "Bad cross product");
#else
    IntVect ivc = cross(iva, ivb);
    if (ivc != IntVect_zero) ++status;
#endif
    IntVect ivd(iva);
    IntVect ive = cross(iva, ivd+1);
    if (ive != IntVect{-1, 2, -1}) ++status;
  }

  // Set (13)
  // Reflect
  {
    static_assert(stc::reflect(IntVect_unit, 0, 0) == IntVect{-1, 1, 1},
                  "Bad reflection");
    static_assert(stc::reflect(IntVect_unit, 1, 0) == IntVect{1, 1, 1},
                  "Bad reflection");
    static_assert(stc::reflect(IntVect_unit, 0, 1) == IntVect{1, -1, 1},
                  "Bad reflection");
    static_assert(stc::reflect(2*IntVect_unit, 0, 2) == IntVect{2, 2, -2},
                  "Bad reflection");
    static_assert(stc::reflect(IntVect{1, 2, 3}, 3, 1) == IntVect{1, 4, 3},
                  "Bad reflection");
    IntVect iva{1,  2,  3};
    iva.reflect(3, 1);
    if (iva != IntVect(1, 4, 3)) ++status;
  }

  // Set (14)
  // Scale
  {
    static_assert(stc::scale(IntVect_unit, 2) == IntVect{ 2, 2, 2 },
                  "Bad scale");
    IntVect iva{ 1, 2, 3 };
    iva.scale(-3);
    if (iva != IntVect{ -3, -6, -9 }) ++status;
  }

  // Set (15)
  // Shift
  {
    static_assert(stc::shift(IntVect_unit, 0, 2) == IntVect{ 3, 1, 1},
                  "Bad shift");
    static_assert(2*stc::shift(IntVect_unit, 0, 2) - 1 == IntVect{ 5, 1, 1},
                  "Bad shift");
    static_assert(stc::shift(IntVect_unit, IntVect{ -1, 1, 2 }) ==
                  IntVect{ 0, 2, 3 }, "Bad shift");
    static_assert(stc::diagShift(IntVect{ -1, 1, 2 }, 1) == IntVect{ 0, 2, 3 },
                  "Bad shift");
    IntVect iva{ 1, 2, 3 };
    iva.shift(0, 2);
    if (iva != IntVect{ 3, 2, 3 }) ++status;
    iva.shift(-2*IntVect_unit);
    if (iva != IntVect{ 1, 0, 1 }) ++status;
    iva.diagShift(1);
    if (iva != IntVect{ 2, 1, 2 }) ++status;
  }

  // Set (16)
  // Magnitude
  {
    if (stc::mag(IntVect_unit) != std::sqrt(c_SpaceDim)) ++status;
    IntVect iva = IntVect_basis(1);
    if (stc::mag(iva) != 1.) ++status;
  }

  // Set (17)
  // Unit vector
  {
    IntVect iva{ 2, 0, 0 };
    if (IntVect(stc::unit(iva)) != IntVect_basis(0)) ++status;
    if (IntVect(stc::unit(3*IntVect_basis(1))) != IntVect_basis(1)) ++status;
  }

  // Set (18)
  // Scalar assignment
  {
    IntVect iva{ 1, 2, 3 };
    iva = -1;
    if (iva != -1*IntVect_unit) ++status;
    iva *= 2;
    if (iva != -2*IntVect_unit) ++status;
    iva += 4;
    if (iva != 2*IntVect_unit) ++status;
    iva /= 2;
    if (iva != IntVect_unit) ++status;
    iva -= 1;
    if (iva != IntVect_zero) ++status;
  }

  // Set (19)
  // Reverse
  {
    IntVect iva;
    iva.reverse({ 3, 2, 1 });
    if (iva != IntVect{ 1, 2, 3 }) ++status;
    iva.reverse();
    if (iva != IntVect{ 3, 2, 1 }) ++status;
    stc::IVec<4> ivb{ 1, 2, 3, 4 };
    if (ivb.reverse() != stc::IVec<4>{ 4, 3, 2, 1 }) ++status;
    if (stc::IVec<4>(ivb).reverse() != stc::IVec<4>{ 1, 2, 3, 4 }) ++status;
    if (ivb != stc::IVec<4>{ 4, 3, 2, 1 }) ++status;
  }

  // Set (20)
  // Stride
  {
    constexpr IntVect iva = stc::cstride(IntVect{ 3, 2, 2 });
    static_assert(iva == IntVect{ 1, 3, 6 }, "Bad stride");
    stc::IVec<4> ivb = stc::stride(stc::IVec<4>{ 3, 2, 2, 2 });
    if (ivb != stc::IVec<4>{ 1, 3, 6, 12 }) ++status;
  }

  // Set (21)
  // Greater than/less than
  {
    static_assert(IntVect{ 2, 1, 1 } < IntVect{ 3, 2, 2 }, "Bad less than");
    static_assert(IntVect{ 3, 2, 2 } > IntVect{ 2, 1, 1 }, "Bad greater than");
    static_assert(IntVect{ 2, 2, 1 } <= IntVect{ 3, 2, 2 },
                  "Bad less than equal");
    static_assert(IntVect{ 3, 2, 2 } >= IntVect{ 2, 2, 1 },
                  "Bad greater than equal");
    static_assert(IntVect_unit <= IntVect_unit, "Bad less than equal");
    static_assert(IntVect_unit >= IntVect_unit, "Bad greater than equal");
  }

  // Set (22)
  // Min/max element direction
  {
    static_assert(stc::indexMinElem(IntVect_unit) == 0, "Bad min index");
    static_assert(stc::indexMinElem(stc::IVec<4>{ 3, 3, 5, 2 }) == 3,
                  "Bad min index");
    static_assert(stc::indexMinElem(stc::IVec<4>{-1, -1, 5, 2 }) == 0,
                  "Bad min index");
    static_assert(stc::indexMinElem(stc::IVec<4>{-1, -2, 5, 2 }) == 1,
                  "Bad min index");
    static_assert(stc::indexMaxElem(IntVect_unit) == 0, "Bad max index");
    static_assert(stc::indexMaxElem(stc::IVec<4>{ 3, 3, -1, 4 }) == 3,
                  "Bad min index");
    static_assert(stc::indexMaxElem(stc::IVec<4>{ 5, 5, -1, 4 }) == 0,
                  "Bad min index");
    static_assert(stc::indexMaxElem(stc::IVec<4>{ 5, 6, -1, 4 }) == 1,
                  "Bad min index");
  }

  // Set (23)
  // Any dim nested loops
  {
    IntVect accum = IntVect_zero;
    IntVect first, second, last;
    const int c0 = 0;
    const int c1 = 1;
    const int cN3 = std::pow(3, c_SpaceDim);
    int c = c0;
    stc::nestedLoop(-IntVect_unit, IntVect_unit,
                    [&]
                    (const IntVect& iv)
                    {
                      accum += iv;
                      if (c == c0)   first = iv;
                      if (c++ == c1) second = iv;
                      if (c == cN3)  last = iv;
                      // if (verbose)
                      //   {
                      //     pout() << iv << std::endl;
                      //   }
                    });
    if (c != cN3) ++status;
    if (first != -IntVect_unit) ++status;
    if (second != IntVect{0, -1, -1}) ++status;
    if (last != IntVect_unit) ++status;
    if (accum != IntVect_zero) ++status;
    if (verbose)
      {
        pout() << "First : " << first << std::endl;
        pout() << "Second: " << second << std::endl;
        pout() << "Last  : " << last << std::endl;
        pout() << "Accum : " << accum << std::endl;
      }
    accum = IntVect_zero;
    c = c0;
    stc::nestedLoop(IntVect_zero, 2*IntVect_unit,
                    [&]
                    (const IntVect& iv)
                    {
                      accum += iv;
                      if (c == c0)   first = iv;
                      if (c++ == c1) second = iv;
                      if (c == cN3)  last = iv;
                      // if (verbose)
                      //   {
                      //     pout() << iv << std::endl;
                      //   }
                    });
    if (c != cN3) ++status;
    if (first != IntVect_zero) ++status;
    if (second != IntVect{1, 0, 0}) ++status;
    if (last != 2*IntVect_unit) ++status;
    if (accum != 27*IntVect_unit) ++status;
    if (verbose)
      {
        pout() << "First : " << first << std::endl;
        pout() << "Second: " << second << std::endl;
        pout() << "Last  : " << last << std::endl;
        pout() << "Accum : " << accum << std::endl;
      }
    accum = IntVect_zero;
    c = c0;
    const int cN2 = std::pow(2, c_SpaceDim);
    stc::nestedLoop(IntVect_zero, 2*IntVect_unit, 2*IntVect_unit,
                    [&]
                    (const IntVect& iv)
                    {
                      accum += iv;
                      if (c == c0)   first = iv;
                      if (c++ == c1) second = iv;
                      if (c == cN2)  last = iv;
                      // if (verbose)
                      //   {
                      //     pout() << iv << std::endl;
                      //   }
                    });
    if (c != cN2) ++status;
    if (first != IntVect_zero) ++status;
    if (second != IntVect{2, 0, 0}) ++status;
    if (last != 2*IntVect_unit) ++status;
    if (accum != 8*IntVect_unit) ++status;
    if (verbose)
      {
        pout() << "First : " << first << std::endl;
        pout() << "Second: " << second << std::endl;
        pout() << "Last  : " << last << std::endl;
        pout() << "Accum : " << accum << std::endl;
      }
  }

  // Set (24)
  // Constructing from a pointer
  {
    int data[3] = { 1, 2, 3 };
    stc::IVec<3> iv(data);
    stc::IVecAlias<3> iva(data);
    stc::IVecConstAlias<3> ivca(data);
    iva[1] = -2;
    if (iv != IntVect{1, 2, 3}) ++status;
    if (iva != IntVect{1, -2, 3}) ++status;
    if (ivca != IntVect{1, -2, 3}) ++status;
    if (verbose)
      {
        pout() << "data: " << data[0] << ' ' << data[1] << ' ' << data[2]
               << std::endl;
        pout() << "iv  : " << iv << std::endl;
        pout() << "iva : " << iva << std::endl;
        pout() << "ivca: " << ivca << std::endl;
      }
    iva = stc::make_IVec<3>::zero();  // Assignment of different types
    if (verbose)
      {
        pout() << "ivca: " << ivca << std::endl;
      }
  }

  // Set (25)
  // Alias new data
  {
    stc::IVec<4> dataA          { -1,  1,  2,  3 };
    constexpr stc::IVec<4> dataB{  4,  5,  6,  7 };
    int dataC[4] =              {  8,  9, 10, 11 };
    const int dataD[4] =        { 12, 13, 14, 15 };
    stc::IVecAlias<3> iv(dataA, 1);
    stc::IVecConstAlias<3> ivc(dataB, 1);
    if (iv  != IntVect{1, 2, 3}) ++status;
    if (ivc != IntVect{5, 6, 7}) ++status;
    if (verbose)
      {
        pout() << "iv : " << iv << std::endl;
        pout() << "ivc: " << ivc << std::endl;
      }
    iv.alias(dataC, 1);
    ivc.alias(dataD, 1);
    if (iv  != IntVect{ 9, 10, 11}) ++status;
    if (ivc != IntVect{13, 14, 15}) ++status;
    if (verbose)
      {
        pout() << "iv : " << iv << std::endl;
        pout() << "ivc: " << ivc << std::endl;
      }
    iv.alias(dataA);
    ivc.alias(dataB);
    if (iv  != IntVect{-1,  1,  2}) ++status;
    if (ivc != IntVect{ 4,  5,  6}) ++status;
    if (verbose)
      {
        pout() << "iv : " << iv << std::endl;
        pout() << "ivc: " << ivc << std::endl;
      }
    dataA[1] = 3;
    if (iv  != IntVect{-1,  3,  2}) ++status;
  }

  // Set (26)
  // Test Morton ordering
  {
    constexpr stc::IVec<2> a{ 0, 0 };
    constexpr stc::IVec<2> b{ 1, 0 };
    constexpr stc::IVec<2> c{ 0, 1 };
    constexpr stc::IVec<2> d{ 1, 1 };
    constexpr stc::IVec<2> e{ 2, 0 };
    constexpr stc::IVec<2> f{ 0, 2 };
    constexpr bool res = stc::morton(b, b);
    if (res) ++status;
    if (!stc::morton(a, b)) ++status;
    if (stc::morton(b, a)) ++status;
    if (!stc::morton(b, c)) ++status;
    if (stc::morton(c, b)) ++status;
    if (!stc::morton(c, d)) ++status;
    if (stc::morton(d, c)) ++status;
    if (!stc::morton(d, e)) ++status;
    if (stc::morton(e, d)) ++status;
    if (!stc::morton(e, f)) ++status;
    if (stc::morton(f, e)) ++status;
    constexpr stc::IVec<3> g{ 1, 1, 0 };
    constexpr stc::IVec<3> h{ 0, 0, 1 };
    constexpr stc::IVec<3> i{ 1, 0, 1 };
    constexpr stc::IVec<3> j{ 0, 1, 1 };
    constexpr stc::IVec<3> k{ 1, 1, 1 };
    if (!stc::morton(g, h)) ++status;
    if (stc::morton(h, g)) ++status;
    if (!stc::morton(h, i)) ++status;
    if (stc::morton(i, h)) ++status;
    if (!stc::morton(i, j)) ++status;
    if (stc::morton(j, i)) ++status;
    if (!stc::morton(j, k)) ++status;
    if (stc::morton(k, j)) ++status;
  }

  // Set (27)
  // Log2 and constructing multi-dim index from linear index
  {
    IntVect ref{ 2, 4, 1 };
    IntVect log2Ref = stc::log2(ref);
    if (log2Ref != IntVect{ 1, 2, 0 }) ++ status;
    for (int i = 0; i != 19; ++i)
      {
        IntVect ivA = stc::fromLinearLog2Dim(i, log2Ref);
        if (ivA != IntVect{ i % 2, i/2 % 4, i/8 }) ++status;
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
