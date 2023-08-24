# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase = turbChem

## location of 'lib' directory
CHOMBO_HOME = ../../..

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames = AMRTools AMRTimeDependent BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ../src

include $(CHOMBO_HOME)/mk/Make.example

# additional libraries
LIBFLAGS += $(LAPACKLIBS) $(FFTWLIBS)
