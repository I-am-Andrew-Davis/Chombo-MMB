# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_ThreadTools

ebase = testHardwareTopology testThreadArchitect \
  threadsEx subTaskTeamsEx 

#ebase = testThreadSupport testHardwareTopology testThreadBuilder       \
#  testThreadArchitect threadsEx subTask

LibNames = ThreadTools BoxTools

include $(CHOMBO_HOME)/mk/Make.test
