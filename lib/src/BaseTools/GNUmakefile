# -*- Mode: Makefile -*- 

### This makefile produces a library for the name in the 'lbase' variable
### depending on the libraries named in the 'LibNames' variable

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

lbase := basetools
cuda_lbase := basetools
LibNames :=

# put this first so 'Cfile' isn't the default target
include $(CHOMBO_HOME)/mk/Make.lib

.PHONY: Cfile
