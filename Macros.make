
# This file is auto-generated, do not edit. If you want to change
# sharedlib flags, you can edit the cmake_macros in this case. You
# can change flags for specific sharedlibs only by checking COMP_NAME.

CFLAGS :=  -qno-opt-dynamic-align -fp-model precise -std=gnu99 -O2 -debug minimal -qopt-report -xCORE_AVX2 -no-fma
CPPDEFS := $(CPPDEFS)  -DCESMCOUPLED -DFORTRANUNDERSCORE -DCPRINTEL
CXX_LDFLAGS :=  -cxxlib
CXX_LINKER := FORTRAN
FC_AUTO_R8 := -r8
FFLAGS :=  -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source -O2 -debug minimal -qopt-report -xCORE_AVX2 -no-fma
FFLAGS_NOOPT := -O0
FIXEDFLAGS := -fixed
FREEFLAGS := -free
HAS_F2008_CONTIGUOUS := TRUE
MACRO_FILE := 
MPICC := mpicc
MPICXX := mpicxx
MPIFC := mpif90
NETCDF_PATH := /glade/u/apps/ch/opt/netcdf-mpi/4.8.1/mpt/2.22/intel/19.1.1/
PIO_FILESYSTEM_HINTS := gpfs
PNETCDF_PATH := /glade/u/apps/ch/opt/pnetcdf/1.12.2/mpt/2.22/intel/19.1.1/
SCC := icc
SCXX := icpc
SFC := ifort
SLIBS := $(SLIBS)  -mkl=cluster
SUPPORTS_CXX := TRUE

ifeq "$(COMP_NAME)" "mom"
  CPPDEFS := $(CPPDEFS)  -DCESMCOUPLED -Duse_LARGEFILE -DFORTRANUNDERSCORE -DCPRINTEL
  FFLAGS :=  $(FC_AUTO_R8)  -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source -O2 -debug minimal -qopt-report -xCORE_AVX2 -no-fma
endif
ifeq "$(COMP_NAME)" "gptl"
  CPPDEFS := $(CPPDEFS)  -DCESMCOUPLED -DFORTRANUNDERSCORE -DCPRINTEL -DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY
endif
