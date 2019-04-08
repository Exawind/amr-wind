
# Parallelisation options
USE_MPI = TRUE
USE_OMP = FALSE

# Debug mode?
DEBUG = FALSE

# Use HYPRE solver?
USE_HYPRE = FALSE

# Profiling
PROFILE       = FALSE
TINY_PROFILE  = FALSE



########################################################################\
# 																		#
# Below are settings which we probably don't want to change very often. #
# 																		#
########################################################################/

# Path to AMReX directory and incflo directories
AMREX_HOME ?= ../../amrex
TOP = ..

# Use OS-friendly compiler
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	COMP = gnu
else ifeq ($(UNAME), Darwin)
	COMP = llvm
endif

# Non-verbose compilation
VERBOSE = FALSE

# Always use 3 dimensions for incflo
DIM = 3

EBASE     ?= incflo

USE_MG        = TRUE
USE_EB        = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

#These are the directories in incflo/src
Bdirs 	:= src
Bdirs 	+= src/boundary_conditions
Bdirs 	+= src/convection
Bdirs 	+= src/derive
Bdirs 	+= src/diffusion
Bdirs 	+= src/embedded_boundaries
Bdirs 	+= src/projection
Bdirs 	+= src/rheology
Bdirs 	+= src/setup
Bdirs 	+= src/utilities

Bpack	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs	+= $(foreach dir, $(Bdirs), $(TOP)/$(dir))

include $(Bpack)
INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

#These are the directories in AMReX
Pdirs   := Base AmrCore Boundary EB
Pdirs   += Extern/Algoim

ifeq ($(USE_HYPRE), TRUE)
Pdirs   += Extern/HYPRE
endif

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
Plocs	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir))

include $(Ppack)
INCLUDE_LOCATIONS += $(Plocs)
VPATH_LOCATIONS   += $(Plocs)

include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/C_CellMG
VPATH_LOCATIONS   += $(AMREX_HOME)/Src/LinearSolvers/C_CellMG

include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/LinearSolvers/MLMG
VPATH_LOCATIONS   += $(AMREX_HOME)/Src/LinearSolvers/MLMG

ALGOIM_HOME ?= $(TOP)/src/extern/algoim
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/Algoim
INCLUDE_LOCATIONS += $(ALGOIM_HOME)/src
VPATH_LOCATIONS   += $(AMREX_HOME)/Src/Extern/Algoim
VPATH_LOCATIONS   += $(ALGOIM_HOME)

INCLUDE_LOCATIONS += $(BLITZ_HOME)/include

DEFINES += -DUSE_ALGOIM

all: $(executable)
	$(SILENT) $(RM) AMReX_buildInfo.cpp
	@echo SUCCESS

# job_info support
CEXE_sources += AMReX_buildInfo.cpp
CEXE_headers += $(AMREX_HOME)/Tools/C_scripts/AMReX_buildInfo.H
INCLUDE_LOCATIONS +=  $(AMREX_HOME)/Tools/C_scripts

AMReX_buildInfo.cpp:
	$(AMREX_HOME)/Tools/C_scripts/makebuildinfo_C.py \
          --amrex_home "$(AMREX_HOME)" \
          --COMP "$(COMP)" --COMP_VERSION "$(COMP_VERSION)" \
          --CXX_comp_name "$(CXX)" --CXX_flags "$(CXXFLAGS) $(CPPFLAGS) $(includes)" \
          --F_comp_name "$(F90)" --F_flags "$(F90FLAGS)" \
          --link_flags "$(LDFLAGS)" --libraries "$(libraries)" \
          --GIT "$(TOP) $(AMREX_HOME)"

vpath %.c   . $(VPATH_LOCATIONS)
vpath %.cpp . $(VPATH_LOCATIONS)
vpath %.h   . $(VPATH_LOCATIONS)
vpath %.H   . $(VPATH_LOCATIONS)
vpath %.F   . $(VPATH_LOCATIONS)
vpath %.f90 . $(VPATH_LOCATIONS)
vpath %.f   . $(VPATH_LOCATIONS)
vpath %.fi  . $(VPATH_LOCATIONS)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) AMReX_buildInfo.cpp

