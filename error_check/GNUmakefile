AMREX_HOME ?= /home/dbeckers/amr-wind/submods/amrex

PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE

#DIM        = 1
#DIM        = 2
DIM        = 3

COMP	   = gcc
USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA = FALSE

BL_NO_FORT = TRUE

EBASE = compareMultilevelToReference

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

CEXE_sources += $(EBASE).cpp

Pdirs 	:= Base Boundary AmrCore Extern/amrdata

Bpack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)
Blocs	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir))

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

CEXE_sources += WritePlotFile.cpp
CEXE_headers += WritePlotFile.cpp
INCLUDE_LOCATIONS += $(AMREX_HOME)/Tools/C_util
VPATH_LOCATIONS += $(AMREX_HOME)/Tools/C_util

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

print-%: ; @echo $* is $($*) ; $(SILENT) $(RM) buildInfo.cpp
