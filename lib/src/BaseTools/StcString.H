
#ifndef _STCSTRING_H_
#define _STCSTRING_H_


/******************************************************************************/
/**
 * \file StcString.H
 *
 * \brief Static string class based on static vector
 *
 *//*+*************************************************************************/


#include <string>

#include "StcVector.H"

#include "BaseNamespaceHeader.H"

namespace stc
{


/******************************************************************************
 */
/// Static string
/**
 *  Primarily used to build fixed-size c-style strings.  Characters added
 *  beyond the fixed-size are silently dropped.  Aliases of existing
 *  stc::String objects can also be build using the VectorAliasImpl and
 *  VectorConstAliasImpl.  There are no constructors or assignments that take
 *  representations so many of the opearators defined for stc::Vector will not
 *  work here.  Operators += and + are redefined as concatenation.
 *
 *  \tparam N           Size of the string (+1 is added for trailing \0).  If
 *                      You want exactly 32 bytes of storage, use N=31
 *  \tparam Impl        Implementation class:
 *                      DefaultImpl          - static storage
 *                      VectorAliasImpl      - alias another stc::String
 *                      VectorConstAliasImpl - const alias another stc::String
 *
 *//*+*************************************************************************/

template <array_size_type N,
          template <typename, array_size_type> class Impl = DefaultImpl>
class String : public VectorImpl<char, N+1, Impl>
{
  using VectorImpl_type = VectorImpl<char, N+1, Impl>;

/*--------------------------------------------------------------------*
 * Types
 *--------------------------------------------------------------------*/

public:

  using is_string = std::true_type;

/*--------------------------------------------------------------------*
 * Constructors and destructors
 *--------------------------------------------------------------------*/

  /// Default constructor
  HOSTDEVICE constexpr String() noexcept
    :
    VectorImpl_type('\0')
    { }

  /// Construct with a single char
  HOSTDEVICE constexpr String(const char a_char) noexcept
    :
    VectorImpl_type('\0')
    {
      static_assert(std::is_same<Impl<char, N>, DefaultImpl<char, N>>::value,
                    "aliasing a char is not permitted");
      if (N > 0)
        {
          VectorImpl_type::m_data[0] = a_char;
        }
    }

  /*
    For making an array of chars that is constexpr.
  */
  struct CharPtr_E
  {
    HOSTDEVICE constexpr CharPtr_E(const char *const a_cstr) noexcept
      :
      m_cstr(a_cstr),
      m_atEnd(false)
      { }
    HOSTDEVICE constexpr char operator()(const array_size_type a_idx)
      noexcept
      {
        if (!m_atEnd && m_cstr[a_idx] == '\0')
          {
            m_atEnd = true;
          }
        return (m_atEnd) ? '\0' : m_cstr[a_idx];
      }
    const char *const m_cstr;
    bool m_atEnd;
  };

  /// Construct from a std::string
  constexpr String(const std::string& a_str) noexcept
    :
    // This is redundant but we need some constructor that is constexpr if
    // we want this object to be constexpr.  VectorImpl_type() cannot be
    // constexpr because DefaultImpl cannot have a constructor defined since
    // it is an aggregate.
    VectorImpl_type('\0')
    {
      static_assert(std::is_same<Impl<char, N>, DefaultImpl<char, N>>::value,
                    "aliasing a std::string is not permitted");
      VectorImpl_type::m_data =
        make_array<char, N+1>(CharPtr_E(a_str.c_str()));
      VectorImpl_type::m_data[N] = '\0';
    }

  /// Construct from char array pointer
  HOSTDEVICE constexpr String(const char *const a_cstr) noexcept
    :
    VectorImpl_type('\0')
    {
      static_assert(std::is_same<Impl<char, N>, DefaultImpl<char, N>>::value,
                    "aliasing a char array is not permitted");
      VectorImpl_type::m_data =
        make_array<char, N+1>(CharPtr_E(a_cstr));
      VectorImpl_type::m_data[N] = '\0';
    }

//--Alias (via constructor): if source has storage and destination is alias

  /// Alias a const vector
  template <array_size_type N2,
            template <typename, array_size_type> class Impl2,
            std::enable_if_t<
              (N2 >= N)
              &&
              std::is_same<Impl<char, N2>,
                           VectorConstAliasImpl<char, N2>>::value, int> = 0
            >
  HOSTDEVICE constexpr String(const String<N2, Impl2>& a_other,
                              const array_size_type    a_idxBeg = 0) noexcept
    :
    VectorImpl_type(a_other, a_idxBeg)
    { 
      assert(a_other.length() - a_idxBeg < N);
    }

  /// Alias a non-const vector
  template <array_size_type N2,
            template <typename, array_size_type> class Impl2,
            std::enable_if_t<
              (N2 >= N)
              &&
              std::is_same<Impl<char, N2>, VectorAliasImpl<char, N2>>::value,
              int> = 0
            >
  HOSTDEVICE String(String<N2, Impl2>&    a_other,
                    const array_size_type a_idxBeg = 0) noexcept
    :
    VectorImpl_type(a_other, a_idxBeg)
    {
      assert(a_other.length() - a_idxBeg < N);
    }

  /// Offset copy
  template <array_size_type N2,
            template <typename, array_size_type> class Impl2,
            std::enable_if_t<
              (N2 >= N)
              &&
              std::is_same<Impl<char, N2>, DefaultImpl<char, N2>>::value,
              int> = 0
            >
  HOSTDEVICE constexpr String(const String<N2, Impl2>& a_other,
                              const array_size_type    a_idxBeg) noexcept
    :
    VectorImpl_type(a_other, a_idxBeg)
    { 
      VectorImpl_type::m_data[N] = '\0';
    }

//--Copy and assignment from string of different size

  /// Construct from a static string with a different size
  template <array_size_type N2,
            template <typename, array_size_type> class Impl2,
            std::enable_if_t<
              std::is_same<Impl<char, N2>, DefaultImpl<char, N2>>::value,
              int> = 0
            >
  HOSTDEVICE constexpr String(const String<N2, Impl2>& a_other) noexcept
    :
    VectorImpl_type('\0')
    {
      VectorImpl_type::m_data =
        make_array<char, N+1>(CharPtr_E(a_other.c_str()));
      VectorImpl_type::m_data[N] = '\0';
    }

  /// Assignment from a static string with a different size
  template <array_size_type N2,
            template <typename, array_size_type> class Impl2,
            std::enable_if_t<
              std::is_same<Impl<char, N2>, DefaultImpl<char, N2>>::value,
              int> = 0
            >
  HOSTDEVICE String& operator=(const String<N2, Impl2>& a_other) noexcept
    {
      VectorImpl_type::m_data =
        make_array<char, N+1>(CharPtr_E(a_other.c_str()));
      VectorImpl_type::m_data[N] = '\0';
      return *this;
    }

//--Note: there are no constructors or assignments that take representations!

//--Standard

  /// Copy constructor
  String(const String&) = default;
  /// Move constructor
  String(String&&) = default;
  /// Assignment constructor
  String& operator=(const String&) = default;
  /// Move assignment constructor
  String& operator=(String&&) = default;
  /// Destructor
  ~String() = default;

/*--------------------------------------------------------------------*
 * Member functions
 *--------------------------------------------------------------------*/

  /// Number of characters (note: different from size())
  HOSTDEVICE array_size_type length() const noexcept
    {
      return std::strlen(c_str());
    }

//--Assignment

  String& operator=(const char a_char) noexcept
    {
      VectorImpl_type::operator=('\0');
      if (N > 0)
        {
          VectorImpl_type::m_data[0] = a_char;
        }
      return *this;
    }

  String& operator=(const std::string& a_str) noexcept
    {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
      std::strncpy(c_str(), a_str.c_str(), N);
#pragma GCC diagnostic pop
      VectorImpl_type::m_data[N] = '\0';
      return *this;
    }

  String& operator=(const char *const& a_cstr) noexcept
    {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
      std::strncpy(c_str(), a_cstr, N);
#pragma GCC diagnostic pop
      VectorImpl_type::m_data[N] = '\0';
      return *this;
    }

//--Concatenation

  template <array_size_type N2,
            template <typename, array_size_type> class Impl2>
  HOSTDEVICE String& operator+=(const String<N2, Impl2>& a_other) noexcept
    {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
      std::strncat(c_str(), a_other.c_str(), N - length());
#pragma GCC diagnostic pop
      return *this;
    }
  
  HOSTDEVICE String& operator+=(const char a_char) noexcept
    {
      const array_size_type pos = length();
      if (pos < N)
        {
          VectorImpl_type::m_data[pos] = a_char;
          VectorImpl_type::m_data[pos+1] = '\0';
        }
      return *this;
    }

  HOSTDEVICE String& operator+=(const std::string& a_str) noexcept
    {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
      std::strncat(c_str(), a_str.c_str(), N - length());
#pragma GCC diagnostic pop
      return *this;
    }

  HOSTDEVICE String& operator+=(const char *const a_cstr) noexcept
    {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
      std::strncat(c_str(), a_cstr, N - length());
#pragma GCC diagnostic pop
      return *this;
    }

//--Access

  /// Constant data pointer
  HOSTDEVICE constexpr const char* c_str() const noexcept
    {
      return VectorImpl_type::data();
    }

  /// Data pointer
  HOSTDEVICE char* c_str() noexcept
    {
      return VectorImpl_type::data();
    }
};


/*******************************************************************************
 *
 * Class String: external operators
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*
 * + is concatentation
 *--------------------------------------------------------------------*/

/*
  Be careful with these because temporaries may be sized smaller than the final
  destination.  So concatentating two strings may not result in the full
  concatentation.
 */

/// Concatenate two stc::Strings
template <array_size_type N1,
          template <typename, array_size_type> class Impl1,
          array_size_type N2,
          template <typename, array_size_type> class Impl2>
HOSTDEVICE auto
operator+(const String<N1, Impl1>& a_x, const String<N2, Impl2>& a_y) noexcept
{
  constexpr array_size_type I = cmax(N1, N2);
  String<I, DefaultImpl> str(a_x);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
  std::strncat(str.c_str(), a_y.c_str(), I - a_x.length());
#pragma GCC diagnostic pop
  return str;
}

/// Concatenate a stc::String and a char
template <array_size_type N,
          template <typename, array_size_type> class Impl>
HOSTDEVICE String<N, DefaultImpl>
operator+(const String<N, Impl>& a_x, const char a_y) noexcept
{
  String<N, DefaultImpl> str(a_x);
  const array_size_type pos = a_x.length();
  if (pos < N)
    {
      str[pos] = a_y;
      str[pos+1] = '\0';
    }
  return str;
}

/// Concatenate a char and a stc::String 
template <array_size_type N,
          template <typename, array_size_type> class Impl>
HOSTDEVICE String<N, DefaultImpl>
operator+(const char a_x, const String<N, Impl>& a_y) noexcept
{
  String<N, DefaultImpl> str(a_x);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
  strncat(str.c_str(), a_y.c_str(), N-1);
#pragma GCC diagnostic pop
  return str;
}

/// Concatenate a stc::String and a std::string
template <array_size_type N,
          template <typename, array_size_type> class Impl>
String<N, DefaultImpl>
operator+(const String<N, Impl>& a_x, const std::string& a_y) noexcept
{
  String<N, DefaultImpl> str(a_x);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
  strncat(str.c_str(), a_y.c_str(), N - a_x.length());
#pragma GCC diagnostic pop
  return str;
}

/// Concatenate a std::string and a stc::String 
template <array_size_type N,
          template <typename, array_size_type> class Impl>
String<N, DefaultImpl>
operator+(const std::string& a_x, const String<N, Impl>& a_y) noexcept
{
  String<N, DefaultImpl> str(a_x);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
  strncat(str.c_str(), a_y.c_str(), N - a_x.length());
#pragma GCC diagnostic pop
  return str;
}

/// Concatenate a stc::String and a c-style string
template <array_size_type N,
          template <typename, array_size_type> class Impl>
HOSTDEVICE String<N, DefaultImpl>
operator+(const String<N, Impl>& a_x, const char *const a_y) noexcept
{
  String<N, DefaultImpl> str(a_x);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
  strncat(str.c_str(), a_y, N - a_x.length());
#pragma GCC diagnostic pop
  return str;
}

/// Concatenate a c-style string and a stc::String 
template <array_size_type N,
          template <typename, array_size_type> class Impl>
HOSTDEVICE String<N, DefaultImpl>
operator+(const char *const a_x, const String<N, Impl>& a_y) noexcept
{
  String<N, DefaultImpl> str(a_x);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-truncation"
  strncat(str.c_str(), a_y.c_str(), N - std::strlen(a_x));
#pragma GCC diagnostic pop
  return str;
}

/*--------------------------------------------------------------------*
 * Equivalence
 *--------------------------------------------------------------------*/

/// Equivalence between two stc::Strings
/** The embedded c-style strings must compare exactly up to and
 *  including the first terminator.
 */
template <array_size_type N1,
          template <typename, array_size_type> class Impl1,
          array_size_type N2,
          template <typename, array_size_type> class Impl2>
HOSTDEVICE constexpr bool
operator==(const String<N1, Impl1>& a_x, const String<N2, Impl2>& a_y) noexcept
{
  constexpr array_size_type I = cmin(N1, N2);
  for (int i = 0; i != I+1; ++i)
    {
      const char chx = a_x[i];
      if (chx != a_y[i]) return false;
      if (chx == '\0') break;
    }
  return true;
}

/// Equivalence between a stc::String and a std::string
/** Only the first N characters of the std::string are compared.
 *  Otherwise, if the length of the std::string is less than N,
 *  the strings must compare exactly.
 */
template <array_size_type N,
          template <typename, array_size_type> class Impl>
constexpr bool
operator==(const String<N, Impl>& a_x, const std::string& a_y) noexcept
{
  for (int i = 0; i != N; ++i)
    {
      const char chx = a_x[i];
      if (chx != a_y.c_str()[i]) return false;
      if (chx == '\0') break;
    }
  return true;
}

/// Equivalence between a std::string and a stc::String
/** Only the first N characters of the std::string are compared.
 *  Otherwise, if the length of the std::string is less than N,
 *  the strings must compare exactly.
 */
template <array_size_type N,
          template <typename, array_size_type> class Impl>
constexpr bool
operator==(const std::string& a_x, const String<N, Impl>& a_y) noexcept
{
  for (int i = 0; i != N; ++i)
    {
      const char chy = a_y[i];
      if (chy != a_x.c_str()[i]) return false;
      if (chy == '\0') break;
    }
  return true;
}

/// Equivalence a stc::String and a c-style string
/** Only the first N characters of the c-style string are compared.
 *  Otherwise, if the length of the c-style string is less than N,
 *  the strings must compare exactly.
 */
template <array_size_type N,
          template <typename, array_size_type> class Impl>
HOSTDEVICE constexpr bool
operator==(const String<N, Impl>& a_x, const char *const a_y) noexcept
{
  for (int i = 0; i != N; ++i)
    {
      const char chx = a_x[i];
      if (chx != a_y[i]) return false;
      if (chx == '\0') break;
    }
  return true;
}

/// Equivalence between a c-style string and a stc::String
/** Only the first N characters of the c-style string are compared.
 *  Otherwise, if the length of the c-style string is less than N,
 *  the strings must compare exactly.
 */
template <array_size_type N,
          template <typename, array_size_type> class Impl>
HOSTDEVICE constexpr bool
operator==(const char *const a_x, const String<N, Impl>& a_y) noexcept
{
  for (int i = 0; i != N; ++i)
    {
      const char chy = a_y[i];
      if (chy != a_x[i]) return false;
      if (chy == '\0') break;
    }
  return true;
}

/*--------------------------------------------------------------------*
 * Non-equivalence
 *--------------------------------------------------------------------*/

/// Non-equivalence between two stc::Strings
/** The embedded c-style strings must compare exactly up to and
 *  including the first terminator for this to return false.
 */
template <array_size_type N1,
          template <typename, array_size_type> class Impl1,
          array_size_type N2,
          template <typename, array_size_type> class Impl2>
HOSTDEVICE constexpr bool
operator!=(const String<N1, Impl1>& a_x, const String<N2, Impl2>& a_y) noexcept
{
  constexpr array_size_type I = cmin(N1, N2);
  for (int i = 0; i != I+1; ++i)
    {
      const char chx = a_x[i];
      if (chx != a_y[i]) return true;
      if (chx == '\0') break;
    }
  return false;
}

/// Non-equivalence between a stc::String and a std::string
/** Only the first N characters of the std::string are compared.
 *  Otherwise, if the length of the std::string is less than N,
 *  the strings must compare exactly.
 */
template <array_size_type N,
          template <typename, array_size_type> class Impl>
constexpr bool
operator!=(const String<N, Impl>& a_x, const std::string& a_y) noexcept
{
  for (int i = 0; i != N; ++i)
    {
      const char chx = a_x[i];
      if (chx != a_y.c_str()[i]) return true;
      if (chx == '\0') break;
    }
  return false;
}

/// Non-equivalence between a std::string and a stc::String
/** Only the first N characters of the std::string are compared.
 *  Otherwise, if the length of the std::string is less than N,
 *  the strings must compare exactly.
 */
template <array_size_type N,
          template <typename, array_size_type> class Impl>
constexpr bool
operator!=(const std::string& a_x, const String<N, Impl>& a_y) noexcept
{
  for (int i = 0; i != N; ++i)
    {
      const char chy = a_y[i];
      if (chy != a_x.c_str()[i]) return true;
      if (chy == '\0') break;
    }
  return false;
}

/// Non-equivalence between a stc::String and a c-style string
/** Only the first N characters of the c-style string are compared.
 *  Otherwise, if the length of the c-style string is less than N,
 *  the strings must compare exactly.
 */
template <array_size_type N,
          template <typename, array_size_type> class Impl>
HOSTDEVICE constexpr bool
operator!=(const String<N, Impl>& a_x, const char *const a_y) noexcept
{
  for (int i = 0; i != N; ++i)
    {
      const char chx = a_x[i];
      if (chx != a_y[i]) return true;
      if (chx == '\0') break;
    }
  return false;
}

/// Non-equivalence between a c-style string and a stc::String
/** Only the first N characters of the c-style string are compared.
 *  Otherwise, if the length of the c-style string is less than N,
 *  the strings must compare exactly.
 */
template <array_size_type N,
          template <typename, array_size_type> class Impl>
HOSTDEVICE constexpr bool
operator!=(const char *const a_x, const String<N, Impl>& a_y) noexcept
{
  for (int i = 0; i != N; ++i)
    {
      const char chy = a_y[i];
      if (chy != a_x[i]) return true;
      if (chy == '\0') break;
    }
  return false;
}


/*******************************************************************************
 *
 * Class String: external related functions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
/// Output the string
/** \param[in] a_os     The output stream
 *  \param[in] a_str    String to output
 *  \return             The output stream
 *//*-----------------------------------------------------------------*/

template <stc::array_size_type N,
          template<typename, stc::array_size_type> class Impl>
inline std::ostream&
operator<<(std::ostream&               a_os,
           const stc::String<N, Impl>& a_str)
{
  a_os.write(a_str.c_str(), a_str.length());
  return a_os;
}

/*--------------------------------------------------------------------*/
//  Input the string
/** \param[in] a_os     The input stream
 *  \param[in] a_str    String to load
 *  \return             The input stream
 *
 *  This behaves identically to reading into a std::string and then
 *  assigning that string to a_str.  If reading to the std::string
 *  fails, a_str is not modified.
 *//*-----------------------------------------------------------------*/

template <stc::array_size_type N,
          template<typename, stc::array_size_type> class Impl>
inline std::istream&
operator>>(std::istream& a_is, stc::String<N, Impl>& a_str)
{
  std::string strtmp;
  a_is >> strtmp;
  if (a_is)
    {
      a_str = strtmp;
    }
  return a_is;
}

}  // namespace stc

#include "BaseNamespaceFooter.H"

#endif  /* ! defined _STCSTRING_H_ */
