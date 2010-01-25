SHELL = /bin/bash # For our sanity!

#-----
# Compiler configuration

# pre-processing.
CPP = cpp -P -C -traditional
CPPFLAGS = -D__Linux -DLAPACK -DFFT_DEFAULT -DPOINTER8 -DMAXMEM='$(MAXMEM)' -D_VCS_VER='$(VCS_VERSION)' $(WORKING_DIR_CHANGES) 

# use compiler with perl scripts to avoid cascade compilation.
compiler = INTEL-ifort9-on-LINUX

# fortran compiler and flags.
FC = gfortran
FFLAGS = -m64 -fcray-pointer
FNEWFLAGS = 
F90FLAGS = -ffree-line-length-none

# flag to tell fortran compiler where to search for module files.
MODULEFLAG = -J

# c compiler and flags.
CC = gcc
CFLAGS = -m64 -O2

# linker, linker flags and libraries.
LD = gfortran
LDFLAGS = -m64
LFLAGS = -L~/local/lib/fftw-3.2/64/gfortran/lib -lfftw3 -framework Accelerate

# For building neci library.
AR = ar
ARFLAGS = -rcs

#-----
# Directory structure and setup

# Directories containing source files (space separated list).
SRC = src

# Directories in which compiled objects are placed.
DEST = dest/opt
# REAL (molecular and gamma-point) objects
GDEST = $(DEST)/real
# COMPLEX (k-point) objects
KDEST = $(DEST)/cmplx

# Directory for compiled executable.
EXE = bin

# Directory for neci libraries.
LIB = lib

# Directory containing scripts to avoid cascade compilation
TOOLS = tools

# Directories that make searches for pre-requisites in addition to the local directory (./).
# Colon separated list.
# We include $(DEST) in this so we can do pre-processing and compilation in 2 steps.
empty :=
space := $(empty) $(empty) # stupid make language...
VPATH = $(subst $(space),:,$(SRC) $(DEST))

# Create output directories if they don't exist.
make_gdest := $(shell test -e $(GDEST) || mkdir -p $(GDEST))
make_kdest := $(shell test -e $(KDEST) || mkdir -p $(KDEST))
make_exe := $(shell test -e $(EXE) || mkdir -p $(EXE))
make_lib := $(shell test -e $(LIB) || mkdir -p $(LIB))

#-----
# VCS info.

# Get the version control id.  Git only.  Outputs a string.
VCS_VERSION := $(shell set -o pipefail && echo -n \" && ( git log --max-count=1 --pretty=format:%H || echo -n 'Not under version control.' ) 2> /dev/null | tr -d '\r\n'  && echo -n \")

# Test to see if the working directory contains changes.  Git only.  If the
# working directory contains changes (or is not under version control) then
# the _WORKING_DIR_CHANGES flag is set.
WORKING_DIR_CHANGES := $(shell git diff --quiet --cached && git diff --quiet 2> /dev/null || echo -n "-D_WORKING_DIR_CHANGES")

#-----
# Additional info.

MAXMEM = 2048 # RAM available, in MB.  Used by the memory logger.

#-----
# Find source files and resultant object files.

# Source extensions.
EXTS = .F90 .F .C .c

# Source filenames.
find_files = $(wildcard $(1)/*$(2))
FSRCFILES := $(foreach dir,$(SRC),$(call find_files,$(dir),.F))
F90SRCFILES := $(foreach dir,$(SRC),$(call find_files,$(dir),.F90))
CSRCFILES := $(foreach dir,$(SRC),$(call find_files,$(dir),.C))
cSRCFILES := $(foreach dir,$(SRC),$(call find_files,$(dir),.c))
SRCFILES := $(FSRCFILES) $(F90SRCFILES) $(CSRCFILES) $(cSRCFILES)

# Objects (strip path and replace extension of source files with .o).
FOBJ_bare := $(addsuffix .o,$(basename $(notdir $(FSRCFILES))))
F90OBJ_bare := $(addsuffix .o,$(basename $(notdir $(F90SRCFILES))))
COBJ_bare := $(addsuffix .o,$(basename $(notdir $(CSRCFILES))))
cOBJ_bare := $(addsuffix .o,$(basename $(notdir $(cSRCFILES))))

# Full path to all objects (gamma-point).
FOBJ := $(addprefix $(GDEST)/, $(FOBJ_bare))
F90OBJ := $(addprefix $(GDEST)/, $(F90OBJ_bare))
COBJ := $(addprefix $(GDEST)/, $(COBJ_bare))
cOBJ := $(addprefix $(GDEST)/, $(cOBJ_bare))
OBJECTS := $(FOBJ) $(F90OBJ) $(COBJ) $(cOBJ)

# Similarly for the kpoint objects.
KFOBJ := $(addprefix $(KDEST)/, $(FOBJ_bare))
KF90OBJ := $(addprefix $(KDEST)/, $(F90OBJ_bare))
KCOBJ := $(addprefix $(KDEST)/, $(COBJ_bare))
KcOBJ := $(addprefix $(KDEST)/, $(cOBJ_bare))

# Objects for standalone NECI.
# We don't need libstub.
OBJECTS_NECI := $(filter-out %libstub.o,$(OBJECTS))

# Objects for CPMD library.
# We don't need necimain.o, cpmdstub.o, init_coul.o, init_could2D.o.  We keep libstub though.
OBJECTS_RCPMD := $(filter-out %necimain.o %cpmdstub.o %init_coul.o %init_could2D.o,$(OBJECTS)) 
OBJECTS_KCPMD := $(addprefix $(KDEST)/,$(notdir $(OBJECTS_RCPMD)))

# Objects for VASP library.
# We don't need necimain.o, vaspstub.o, init_coul.o, init_could2D.o.  We keep libstub though.
OBJECTS_RVASP := $(filter-out %necimain.o %vaspstub.o % %init_coul.o %init_could2D.o, $(OBJECTS)) 
OBJECTS_KVASP := $(addprefix $(KDEST)/,$(notdir $(OBJECTS_RVASP)))

#-----
# Dependency files.

# Fortran dependencies.
# We need these before compiling any fortran files.
FRDEPEND = $(GDEST)/.depend
FCDEPEND = $(KDEST)/.depend
FDEPEND = $(FRDEPEND) $(FCDEPEND)

# C dependencies.
# We don't need these when we first compile, only when we recompile.
# We achieve this (most of the time) by recompiling the C dependencies every
# time we compile.
CDEPEND_FILES = $(COBJ:.o=.d)
cDEPEND_FILES = $(cOBJ:.o=.d)
KCDEPEND_FILES = $(KCOBJ:.o=.d)
KcDEPEND_FILES = $(KcOBJ:.o=.d)
CDEPEND = $(CDEPEND_FILES) $(cDEPEND_FILES) $(KCDEPEND_FILES) $(KcDEPEND_FILES)

#-----
# Goals

.PHONY: clean depend help neci.x

# First, some helpful macros.

# Rebuild of environment_report.
GBLD_ENV = rm $(GDEST)/environment_report.* && $(MAKE) $(GDEST)/environment_report.o
KBLD_ENV = rm $(KDEST)/environment_report.* && $(MAKE) $(KDEST)/environment_report.o

# Creating an archive from *.o files.
ARCHIVE = $(AR) $(ARFLAGS) $@ $^

$(EXE)/neci.x : $(OBJECTS_NECI)
	$(GBLD_ENV)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS_NECI) $(LFLAGS)

neci.x: $(EXE)/neci.x

$(LIB)/gneci-cpmd.a : $(OBJECTS_RCPMD)
	$(GBLD_ENV)
	$(ARCHIVE)

$(LIB)/kneci-cpmd.a : $(OBJECTS_KCPMD)
	$(KBLD_ENV)
	$(ARCHIVE)

$(LIB)/gneci-vasp.a : $(OBJECTS_RVASP)
	$(GBLD_ENV)
	$(ARCHIVE)

$(LIB)/kneci-vasp.a : $(OBJECTS_KVASP)
	$(KBLD_ENV)
	$(ARCHIVE)

clean: 
	  rm -f {$(GDEST),$(KDEST)}/*.{f,f90,mod,o,c,x,a,d} $(EXE)/neci.x $(LIB)/*.a

# Generate dependency files.
$(FRDEPEND):
	$(TOOLS)/sfmakedepend --file - --silent $(SRCFILES) --objdir \$$\(GDEST\) --moddir \$$\(GDEST\) > $(FRDEPEND)

$(FCDEPEND):
	$(TOOLS)/sfmakedepend --file - --silent $(SRCFILES) --objdir \$$\(KDEST\) --moddir \$$\(KDEST\) > $(FCDEPEND)

depend: 
	$(MAKE) -B $(FDEPEND) $(CDEPEND)

#-----
# Compilation macros (explicit rules)

.SUFFIXES:
.SUFFIXES: $(EXTS) .f .f90

# Some more helpful macros.
CPP_BODY = $(CPPFLAGS) $< $@
C_BODY = $(CFLAGS) -c $< -o $@ 
MAKE_C_GDEPS = $(CC) $(CFLAGS) -MM -MT \$$\(GDEST\)/$(addsuffix .o,$(basename $(notdir $@))) $< -o $@
MAKE_C_KDEPS = $(CC) $(CFLAGS) -MM -MT \$$\(KDEST\)/$(addsuffix .o,$(basename $(notdir $@))) $< -o $@

# Compiling fortran source files...

# 1. Pre-processing.
# a) gamma-point.
$(GDEST)/%.f90: %.F90
	$(CPP) $(CPP_BODY)

$(GDEST)/%.f: %.F
	$(CPP) $(CPP_BODY)

# b) k-point.
$(KDEST)/%.f90: %.F90
	$(CPP) -D__CMPLX $(CPP_BODY)

$(KDEST)/%.f: %.F
	$(CPP) -D__CMPLX $(CPP_BODY)

# 2. Compile.
$(F90OBJ) $(KF90OBJ): %.o: %.f90
	perl -w $(TOOLS)/compile_mod.pl -cmp "perl -w $(TOOLS)/compare_module_file.pl -compiler $(compiler)" -fc "$(FC) $(FFLAGS) $(F90FLAGS) $(MODULEFLAG) $(dir $@) -I $(SRC) -c $< -o $@" -provides "$@" -requires "$^"

$(FOBJ) $(KFOBJ): %.o: %.f
	perl -w $(TOOLS)/compile_mod.pl -cmp "perl -w $(TOOLS)/compare_module_file.pl -compiler $(compiler)" -fc "$(FC) $(FFLAGS) $(FNEWFLAGS) $(MODULEFLAG) $(dir $@) -I $(SRC) -c $< -o $@" -provides "$@" -requires "$^"

# Compiling C source files...
# a) gamma-point.
$(COBJ): $(GDEST)/%.o: %.C
	$(CC) $(CPPFLAGS) $(C_BODY)

$(cOBJ): $(GDEST)/%.o: %.c
	$(CC) $(C_BODY)

# b) k-point.
$(KCOBJ): $(KDEST)/%.o: %.C
	$(CC) $(CPPFLAGS) -D__CMPLX $(C_BODY)

$(KcOBJ): $(KDEST)/%.o: %.c
	$(CC) $(C_BODY)

# Update C dependency files.
# a) gamma-point.
$(cDEPEND_FILES): $(GDEST)/%.d: %.c
	$(MAKE_C_GDEPS)

$(CDEPEND_FILES): $(GDEST)/%.d: %.C
	$(MAKE_C_GDEPS)

# b) k-point.
$(KcDEPEND_FILES): $(KDEST)/%.d: %.c
	$(MAKE_C_KDEPS)

$(KCDEPEND_FILES): $(KDEST)/%.d: %.C
	$(MAKE_C_KDEPS)

#-----
# Include dependency files

# Dependency files will be generated if they don't exist.
include $(FDEPEND)
include $(CDEPEND)
