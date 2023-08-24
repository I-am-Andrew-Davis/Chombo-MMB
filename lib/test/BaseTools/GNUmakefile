# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_BaseTools

ebase =  clock testTask testCH_Attach testRefCountedPtr \
   testRefCountedPtrConstruct testParmParse test_complex testRootSolver \
   testArrayAllocator testStcVector testStcString testDynArray testHash \
   testComputeCache testNanoflann testHashTable testCudaVectorAdd \
   testGlobalStruct testCudaSharedPtr

# note that BaseTools library should be included by default, even 
# if we don't specify it here
LibNames =
LibCudaNames = BaseTools

include $(CHOMBO_HOME)/mk/Make.test
