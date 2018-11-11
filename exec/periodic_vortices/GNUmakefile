# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../amrex

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# Use OS-friendly compiler
UNAME  := $(shell uname)
ifeq ($(UNAME), Linux)
	COMP = gnu
else ifeq ($(UNAME), Darwin)
	COMP = llvm
endif

# Parallelisation options
USE_MPI = TRUE
USE_OMP = FALSE

# Debug mode (much slower)
DEBUG     = TRUE

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = FALSE

# physics
DIM      = 3

include ../Make.incflo
