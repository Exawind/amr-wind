# Path to AMReX directory
AMREX_HOME ?= ../../../amrex

# Path to incflo directory
TOP = ../..

# Use OS-friendly compiler
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	COMP = gnu
else ifeq ($(UNAME), Darwin)
	COMP = llvm
endif

# Parallelisation options
USE_MPI = TRUE
USE_OMP = FALSE

# Debug mode (much slower)
DEBUG = TRUE

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = FALSE

# Physics
DIM      = 3

include ../Make.incflo
