
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
HYPRE_DIR ?= ../../hypre/src/hypre

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
Pdirs   := Base AmrCore Boundary EB LinearSolvers/MLMG LinearSolvers/Projections

ifeq ($(USE_HYPRE), TRUE)
Pdirs   += Extern/HYPRE
endif

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

all: $(executable)
	$(SILENT) $(RM) AMReX_buildInfo.cpp
	@echo SUCCESS

# job_info support
CEXE_sources += AMReX_buildInfo.cpp
CEXE_headers += AMReX_buildInfo.H
VPATH_LOCATIONS   +=  $(AMREX_HOME)/Tools/C_scripts
INCLUDE_LOCATIONS +=  $(AMREX_HOME)/Tools/C_scripts

AMReX_buildInfo.cpp:
	$(AMREX_HOME)/Tools/C_scripts/makebuildinfo_C.py \
          --amrex_home "$(AMREX_HOME)" \
          --COMP "$(COMP)" --COMP_VERSION "$(COMP_VERSION)" \
          --CXX_comp_name "$(CXX)" --CXX_flags "$(CXXFLAGS) $(CPPFLAGS) $(includes)" \
          --F_comp_name "$(F90)" --F_flags "$(F90FLAGS)" \
          --link_flags "$(LDFLAGS)" --libraries "$(libraries)" \
          --GIT "$(TOP) $(AMREX_HOME)"

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) AMReX_buildInfo.cpp
