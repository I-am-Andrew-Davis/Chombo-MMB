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
 * \file
 *
 * \brief System dependent functions
 *
 *//*+*************************************************************************/

#include "CH_config.H"

#ifdef CHDEF_SYSTEM_HAVE_POSIXMEMALIGN
#define _XOPEN_SOURCE 600
#endif
#include <unistd.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "CH_System.H"
#include "BaseNamespaceHeader.H"

/*--------------------------------------------------------------------*/
//  Check if a file exists
/**
 *  \param[in]  a_filename
 *                      Name of the file
 *  \return             1 - File exists
 *                      0 - File does not exist
 *//*-----------------------------------------------------------------*/

int CH_System::fileExists(const char *const a_filename)
{
  struct stat buf;
  if (stat(a_filename, &buf) == 0)
    {
      return 1;
    }
  return 0;
}

/*--------------------------------------------------------------------*/
//  Allocate aligned memory
/**
 *  \param[out] a_memptr
 *                      Pointer to allocated memory
 *  \param[in]  a_alignment
 *                      Alignment in bytes.  Must be a multiple of
 *                      sizeof(void*) and a power of 2.
 *  \param[in]  a_size  Number of bytes to allocate
 *  \return             0       - Success
 *                      <posix_memalign>
 *                      EINVAL  - Invalid alignment parameter
 *                      ENOMEM  - Out of memory
 *                      <malloc>
 *                      1       - Out of memory
 *  \note
 *  <ul>
 *    <li> This function returns raw memory.  Use placement new for
 *         construction of objects if required.
 *    <li> Memory allocated with memalign should be deallocated with
 *         free()
 *  </ul>
 *//*-----------------------------------------------------------------*/

int CH_System::memalign(void **a_memptr, size_t a_alignment, size_t a_size)
{
#ifdef CHDEF_SYSTEM_HAVE_POSIXMEMALIGN
    return posix_memalign(a_memptr, a_alignment, a_size);
#else
    *a_memptr = malloc(a_size);
    return (*a_memptr == 0);
#endif
}

/*--------------------------------------------------------------------*/
//  Get the path and name of the currently running executable
/** \param[in]  a_procPath
 *                      Char buffer of length a_len
 *  \param[out] a_procPath
 *                      Path and name of current executable in a C-
 *                      string
 *  \param[in]  a_len   Size of buffer a_procPath
 *  \return             > 0     - Successful -- length of string
 *                                a_procPath
 *                      =-1     - Error
 *                      =-a_len - Likely ran out of space in
 *                                a_procPath
 *//*-----------------------------------------------------------------*/

int CH_System::getProcessPath(char *const a_procPath, const int a_len)
{
  int len = readlink("/proc/self/exe", a_procPath, a_len-1);
  if (len > 0)
    {
      a_procPath[len]='\0';
      if (len == a_len-1)
        {
          len = -a_len;
        }
    }
  return len;
}

#include "BaseNamespaceFooter.H"
