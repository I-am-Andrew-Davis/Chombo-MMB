### This is a compiler-specific definitions file for GNU compilers.

## It sets values for the defaults variables ($def*), which will be used
## if the calling makefile doesn't set the corresponding variables.
## The corresponding variables are the ones without the "def" prefix.

## Defaults variables:
##  defcppcomflags = C-preprocessor options for both C++ and Fortran code
##  defcppdbgflags = C-preprocessor options for both C++ and Fortran code when DEBUG=TRUE
##  defcppoptflags = C-preprocessor options for both C++ and Fortran code when OPT=TRUE or HIGH
##  defcxxcppflags = C-preprocessor options for only C++ code
##  defcxxcomflags = C++ compiler/linker options
##  defcxxdbgflags : for C++ compiling and linking when DEBUG=TRUE
##  defcxxoptflags : for C++ compiling and linking when OPT=TRUE or HIGH
##  defcxxprofflags: for C++ compiling and linking when PROFILE=TRUE
##  deffcppflags   = C-preprocessor options for only Fortran code
##  deffcomflags   = Fortran compiler options
##  deffdbgflag    : for Fortran compiling when DEBUG=TRUE
##  deffoptflag    : for Fortran compiling when OPT=TRUE or HIGH
##  deffprofflags  : for Fortran compiling when PROFILE=TRUE
##  defflibflags   : library options needed to link the Fortran code
##  defldcomflags  = linker options
##  deflddbgflags  = linker options when DEBUG=TRUE
##  defldoptflags  = linker options when OPT=TRUE or HIGH
##  defldprofflags = linker options when PROFILE=TRUE
##  cppcallsfort   = preprocessor macro as used in FORT_PROTO.H to specify
##                    how the C++ code should construct names for calling Fortran
##
## Special object models:
##  USE_64 = TRUE   Compile/link with 64bit pointers and longs
##
##  Note: 64bit only works on AMD64 or PowerPC computers with gcc 3.3 or later.  <dbs,Mar04>
##
## Notes:
##  GNU C/C++ compilers set the preprocessor macros: (use -E -dM options to display)
##    __GNUC__       (int) defined for gcc or g++, == major version number
##    __GNUG__       (int) defined for g++, == major version number
##    __GNUC_MINOR__ (int) defined for gcc or g++, == minor version number

makefiles+=compiler/Make.defs.GNU

################################################################

_gppversion := $(subst ., ,$(shell $(CXX) -dumpversion))
_gppmajorver:= $(word 1,$(_gppversion))
_gppminorver:= $(word 2,$(_gppversion))

# As of at least g++ version 7, -dumpversion only gives the
# major version and -dumpfullversion should be used to get the
# minor version.
ifeq (00,$(shell test $(_gppmajorver) -ge 7; echo $$?)$(shell test -z $(_gppminorver); echo $$?))
  _gppversion := $(subst ., ,$(shell $(CXX) -dumpfullversion))
  _gppminorver:= $(word 2,$(_gppversion))
endif

  # `-ftemplate-depth' is needed to keep mangled function names from
  #  getting too long for some linkers to handle.
  # HDF5 uses `long long', so disable warnings about it.
  # `-Wno-sign-compare' turns off warnings about comparing signed and unsigned int -
  #   this is a meaningful warning, but it appears way too often in our code.
  # `-Wno-unused-local-typedefs' turns off a set of annoying "typedef locally 
  # defined but not used warnings in RefCountedPtr.H
  # enum-compare is noisy in HashIVS

ifeq ($(CXXSTD),11)
  _cxxbaseflags := -std=c++11
else ifeq ($(CXXSTD),14)
  _cxxbaseflags := -std=c++14
else ifeq ($(CXXSTD),17)
  _cxxbaseflags := -std=c++17
else
  _cxxbaseflags := -std=c++11
endif
_cxxbaseflags += -march=native -Wno-long-long -Wno-sign-compare -Wno-deprecated -ftemplate-depth-99 -Wno-unused-local-typedefs -Wno-enum-compare

# Too much noise from this so only enable for debug
ifeq ($(DEBUG),FALSE)
  _cxxbaseflags += -Wno-unused-but-set-variable -Wno-literal-suffix
endif

  # g++ v3.4.0 started objecting to some uses of HOFFSET()
  # TJL - Changed moving from csh to sh...
ifeq (00,$(shell test $(_gppmajorver) -ge 3 ; echo $$?)$(shell test $(_gppminorver) -ge 4 ; echo $$?))
  _cxxbaseflags += -Wno-invalid-offsetof
endif

ifeq (0,$(shell test $(_gppmajorver) -gt 3 ; echo $$?))
  _cxxbaseflags += -Wno-variadic-macros
endif

  # NOTE: -malign-double causes MPI programs that use pout() to
  # segfault in ~ostream(), so dont use it

ifeq ($(OPENMPCC),TRUE)
  _cxxbaseflags += -fopenmp
else
  _cxxbaseflags += -Wno-unknown-pragmas
endif

defcxxcomflags := $(_cxxbaseflags)
defcxxoptflags := -O3
defcxxdbgflags := -g -pedantic -Wall

  ##  # together so force the preprocessor to use the Fortran compiler
ifeq ($(system),Darwin)
  CH_CPP = /usr/bin/cpp -E -P
endif

# the preprocessor in v3 or later strips "//" by default,
# but this is a valid Fortran operator, so force it to leave comments
ifeq (0,$(shell test $(_gppmajorver) -ge 3 ; echo $$?))
  CH_CPP += 
endif

# the preprocessor for v4.8+ adds some header gibberish that we 
# suppress with -nostdinc.
ifeq (00,$(shell test $(_gppmajorver) -ge 4 ; echo $$?)$(shell [ $(_gppminorver) -ge 8 -o $(_gppmajorver) -ge 5 ] ; echo $$?))
  CH_CPP += -nostdinc
endif

################################################################

_g77version := $(subst ., ,$(shell $(FC) -dumpversion))
_g77majorver:= $(word 1,$(_g77version))
_g77minorver:= $(word 2,$(_g77version))

# As of at least gcc 7, -dumpversion only gives the
# major version and -dumpfullversion should be used to get the
# minor version.
ifeq (00,$(shell test $(_gppmajorver) -ge 7; echo $$?)$(shell test -z $(_gfortminorver); echo $$?))
    _gfortversion := $(subst ., ,$(shell $(FC) -dumpfullversion))
    _gfortminorver:= $(word 2,$(_gfortversion))
  endif

deffcomflags = -Wno-line-truncation -fno-second-underscore
# these flags get rid of unused variables warnings. (-Wno-unused-parameter -Wno-unused-variable)
deffcomflags += -Wno-unused-parameter -Wno-unused-variable
deffoptflags = -march=native -O3 -funroll-loops
deffdbgflags = -g -ggdb -Wall -W

  # Try this for gcc 4.2 and above...
  # newer versions of gcc (4.2.1 works) can use this
ifeq (00,$(shell test $(_gppmajorver) -ge 4 ; echo $$?)$(shell test $(_gppminorver) -ge 2 ; echo $$?))
   deffdbgflags += -fbounds-check
endif

  # Note this is a copy/paste from the g77 ifeq above
  # non-default object models are meaningful ONLY on some machines
ifeq ($(USE_64),TRUE)
   # works on an AMD64 machine that defaults to 32bit
   # may also work on a Power{PC,4,5} but I'm not positive
  deffcomflags += -m64
endif

ifeq ($(OPENMPCC),TRUE)
  deffcomflags  += -fopenmp
  defflibflags  += -lgfortran -lm -lgomp
else
  defflibflags  += -lgfortran -lm 
endif

