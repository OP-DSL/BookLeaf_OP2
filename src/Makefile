#
# !Crown Copyright 2014 AWE.
# !
# ! This file is part of Bookleaf.
# !
# ! Bookleaf is free software: you can redistribute it and/or modify it under
# ! the terms of the GNU General Public License as published by the
# ! Free Software Foundation, either version 3 of the License, or (at your option)
# ! any later version.
# !
# ! Bookleaf is distributed in the hope that it will be useful, but
# ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# ! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# ! details.
# !
# ! You should have received a copy of the GNU General Public License along with
# ! Bookleaf. If not, see http://www.gnu.org/licenses/.
#
#
# machine input
MKFILEM = makefile.GENERIC
MKFILEC = makefile.gcc
MKFILEV = $(SRCDIR)/makefiles/makefile.VER
# paths
SRCDIR  = .
ALEDIR  = $(SRCDIR)/ale
COMDIR  = $(SRCDIR)/comms
CUTDIR  = $(SRCDIR)/cutils
EOSDIR  = $(SRCDIR)/eos
HYDDIR  = $(SRCDIR)/hydro
IAODIR  = $(SRCDIR)/io
MODDIR  = $(SRCDIR)/mods
SCPDIR  = $(SRCDIR)/putils
UTLDIR  = $(SRCDIR)/utils

# build info
NO_MPI  =
NO_OMP  =
PHI			=
SILO    =
DEP     =
# Version
VERS = 1.1
EXE  = bookleaf
ifeq ($(MOD),)
  EXENAME  = $(EXE)
else
  EXENAME= $(EXE)_$(MOD)
endif
ifdef OPTARGET
ifeq ($(MOD),)
  EXENAME  = $(EXE)_$(OPTARGET)
else
  EXENAME= $(EXE)_$(MOD)_$(OPTARGET)
endif
endif

# Determine level
ifeq ($(LEV),)
  LEV=OPT
endif
FIND_LEV=$(findstring $(LEV),MEMDEBUG DEBUG TV MEMTV OPTDEBUG RELEASE OPT)
ifeq ($(FIND_LEV),)
  $(error Unknown build level)
endif

# Default
help:
	# Usage
	# Pre_pro:
	#  MKFILEM     - machine specific makefile,  DEFAULT makefile.GENERIC
	#  MKFILEC     - compiler specific makefile, DEFAULT makefile.gcc
	#  MKFILEV     - override default list of library versions
	#  SRCDIR      - path to source directory,   DEFAULT .
	#  LEV         - build level, options: MEMDEBUG,DEBUG,TV,MEMTV,OPTDEBUG,RELEASE,
	#                DEFAULT OPT
	#  MOD         - include problem specific initialisation modification, DEFAULT none
	#  NO_MPI      - Do not compile MPI routines
	#  NO_OMP      - Do not use OpenMP threads
	#  PHI         - Build executable for Intel Xeon Phi
	#  SILO        - Compile Silo routines
	#  TIO         - Compile TyphonIO routines
	#
	#  OPTARGET    - target for OP2 build, either 'genseq', 'openmp', or 'cuda'
	#
	# Commands:
	#  help        - print this message
	#  bookleaf    - build executable
	#  clean       - remove objects, modules and executable

# Include machine specific options
ifneq (,$(findstring $(EXE),$(MAKECMDGOALS)))
  include $(SRCDIR)/makefiles/$(MKFILEM)
  include $(SRCDIR)/makefiles/$(MKFILEC)
  include $(MKFILEV)
endif

# Include object list
ifneq (,$(findstring $(MAKECMDGOALS),$(EXE)))
  include $(SRCDIR)/makefile.OBJ
  include $(ALEDIR)/makefile.ALE
  include $(EOSDIR)/makefile.EOS
  include $(HYDDIR)/makefile.HYDRO
  include $(IAODIR)/makefile.IO
  include $(UTLDIR)/makefile.UTILS
endif

# set pre-processor flags, adjust object list and dependency list
ifneq ($(SILO),)
  ifeq ($(MKFILEC), makefile.xl)
    PP_SILO=-WF,-DSILO
  else
    PP_SILO=-DSILO
  endif
  OBJ+=silo.o
  ifneq (,$(findstring $(MAKECMDGOALS),$(EXE)))
    include $(CUTDIR)/makefile.cutils
  endif
endif
ifneq ($(TIO),)
  ifeq ($(MKFILEC), makefile.xl)
    PP_TIO=-WF,-DTIO
  else
    PP_TIO=-DTIO
  endif
  OBJ+=tio.o
endif
ifneq ($(MOD),)
  ifeq ($(MKFILEC), makefile.xl)
    PP_MOD=-WF,-DMODY
  else
    PP_MOD=-DMODY
  endif
  ifndef OPTARGET
  OBJ+=$(MOD).o
  else
	ifneq ($(OPTARGET),meshgen)
  OBJ+=$(MOD)_op.o
	endif
  endif
	ifneq ($(OPTARGET),meshgen)
  OBJ+=$(MOD)_kernels.o
	endif
endif
ifeq ($(NO_MPI),)
  ifneq (,$(findstring $(MAKECMDGOALS),$(EXE)))
#    include $(COMDIR)/makefile.T3
  endif
else
  ifeq ($(MKFILEC), makefile.xl)
    PP_MPI=-WF,-DNOMPI
  else
    PP_MPI=-DNOMPI
  endif
  ifneq (,$(findstring $(MAKECMDGOALS),$(EXE)))
#    include $(COMDIR)/makefile.T1
  endif
endif
ifneq ($(NO_OMP),)
  ifeq ($(MKFILEC), makefile.xl)
    PP_OMP=-WF,-DNOOMP
  else
    PP_OMP=-DNOOMP
  endif
endif

# Include dependencies
ifneq (,$(findstring $(EXE),$(MAKECMDGOALS)))
  sinclude ./makefile.dep
endif

# Source paths
VPATH=$(ALEDIR):$(CUTDIR):$(EOSDIR):$(HYDDIR):$(IAODIR):      \
      $(MODDIR):$(SRCDIR):$(UTLDIR) #$(COMDIR)
# Collect all pre-processor options
PREPRO += $(PP_MOD)   \
	 $(PP_MPI)   \
	 $(PP_SILO)  \
	 $(PP_TIO)  \
	 $(PP_OMP)

CPREPRO = $(F2C_DEF)

# Includes
INC =	$(SILOINC)   \
	$(HDF5INC)	\
	$(TIOINC)	\
	-I$(OP2_INSTALL_PATH)/fortran/mod/$(OP2_COMPILER)

# Libraries
LIB = 	$(SILOLIB)   \
	$(HDF5LIB)   \
	$(TIOLIB)   \
	$(ZLIBLIB)

# Rules
.SUFFIXES:
.SUFFIXES: .o .f90 .F90 .CUF .c
.f90.o:
	$(FC) $(FFLAGS) $(CUDAF) $(INC) $<
.F90.o:
	$(FC) $(FFLAGS) $(CUDAF) $(PREPRO) $(INC) $<
.c.o:
	$(CC) $(CFLAGS) $(CUDAF) $(CPREPRO) $<
.CUF.o:
	$(FC) $(FFLAGS) $(CUDAF) $(INC) $<
# Dependencies
$(OBJ): | makefile.dep
DEP+=makefile.dep

# Commands
$(EXE): $(DEP) $(OBJ)
	$(FC) $(OBJ) $(LIB) $(LFLAGS) $(EXENAME) -L$(OP2_INSTALL_PATH)/fortran/lib/ $(OP2_LIB)

clean:
	- rm -f *.o
	- rm -f *.oo
	- rm -f *.mod
	- rm -f *.lst
	- rm -f *.ipa
	- rm -f bookleaf*
	- rm -f makefile.dep

makefile.dep:
	$(SCPDIR)/depend.py -d $(SRCDIR) $(OBJ)
