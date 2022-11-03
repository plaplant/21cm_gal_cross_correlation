# HDF5
#H5HOME = /opt/homebrew
H5HOME = /opt/local/hdf5-1.12.1
H5FC   = $(H5HOME)/bin/h5fc
H5INC  = -I$(H5HOME)/include
H5LIB  = -L$(H5HOME)/lib

INC = $(H5INC)
LIB = $(H5LIB)

# compiler
FC = $(H5FC)
# test for gfortran or ifort
compiler_type = $(shell $(FC) --version 2>&1 | head -1 | cut -d " " -f 1)
ifeq ($(compiler_type),ifort)
# ifort options
FFLAGS = -qopenmp -O3 -ip -ipo -shared-intel -fpp -DHAVE_IFORT
DB0    = -g -traceback
DB1    = -check -debug -warn
else
# gfortran options
FFLAGS = -fopenmp -O3 -cpp
DB0    = -g -fbacktrace
DB1    = -ggdb -Wextra -fbounds-check
endif

# object files
OBJ = cosmo_parameters.o physical_constants.o general_tools.o cosmo_tools.o

# debug?
ifeq ($(DEBUG), 1)
FFLAGS += $(DB0) $(DB1)
else ifeq ($(DEBUG), 0)
FFLAGS +=
else
FFLAGS += $(DB0)
endif

# executables
compute_snr.x: $(OBJ) compute_snr.f90
	$(FC) $(FFLAGS) $(INC) $(OBJ) compute_snr.f90 -o $@ $(LIB)

# object files
%.o: %.f90
	$(FC) $(FFLAGS) $(INC) -c $*.f90

.PHONY: all
all: compute_snr.x

.PHONY: clean
clean:
	rm -rf *.o *.mod *.x
