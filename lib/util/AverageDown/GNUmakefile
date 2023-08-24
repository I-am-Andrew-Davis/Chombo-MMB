# -*- Mode: Makefile -*- 

### This makefile produces an executable for each 
### name in the `ebase' variable
ebase = averagedown averagedownModified averageUpModified

## location of 'lib' directory
CHOMBO_HOME = ../..

##
## names of Chombo libraries needed by this program, in order of search.
##
LibNames = HOMappedConservation \
           HOMappedGeometry \
           MOLAMRTimeDependent \
           AMRTimeDependent \
           HOAMRTools \
           AMRTools \
           BoxTools \
           BaseTools


# input file for 'run' target
INPUT = inputs.averagedown
override VERBOSE = # #program doesnt handle -q option


include $(CHOMBO_HOME)/mk/Make.test

LIBFLAGS += $(LAPACKLIBS) $(FFTWLIBS)
