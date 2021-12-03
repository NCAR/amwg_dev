SUPPORTS_CXX := FALSE
ifeq ($(COMPILER),gnu)
  FFLAGS :=   -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none 
  SUPPORTS_CXX := TRUE
  CFLAGS :=  -std=gnu99 
  CXX_LINKER := FORTRAN
  FC_AUTO_R8 :=  -fdefault-real-8 -fdefault-double-8
  FFLAGS_NOOPT :=  -O0 
  FIXEDFLAGS :=   -ffixed-form 
  FREEFLAGS :=  -ffree-form 
  HAS_F2008_CONTIGUOUS := FALSE
  MPICC :=  mpicc  
  MPICXX :=  mpicxx 
  MPIFC :=  mpif90 
  SCC :=  gcc 
  SCXX :=  g++ 
  SFC :=  gfortran 
endif
ifeq ($(COMPILER),intel)
  FFLAGS :=  -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source  
  SUPPORTS_CXX := TRUE
  CFLAGS :=   -qno-opt-dynamic-align -fp-model precise -std=gnu99 
  CXX_LINKER := FORTRAN
  FC_AUTO_R8 :=  -r8 
  FFLAGS_NOOPT :=  -O0 
  FIXEDFLAGS :=  -fixed  
  FREEFLAGS :=  -free 
  MPICC :=  mpicc  
  MPICXX :=  mpicxx 
  MPIFC :=  mpif90 
  SCC :=  icc 
  SCXX :=  icpc 
  SFC :=  ifort 
  CXX_LDFLAGS :=  -cxxlib 
endif
ifeq ($(COMPILER),pgi)
  FFLAGS :=   -i4 -gopt  -time -Mextend -byteswapio -Mflushz -Kieee  
  CFLAGS :=  -gopt  -time 
  CXX_LINKER := CXX
  FC_AUTO_R8 :=  -r8 
  FFLAGS_NOOPT :=  -O0 
  FIXEDFLAGS :=  -Mfixed 
  FREEFLAGS :=  -Mfree 
  HAS_F2008_CONTIGUOUS := FALSE
  LDFLAGS :=  -time -Wl,--allow-multiple-definition 
  MPICC :=  mpicc 
  MPICXX :=  mpicxx 
  MPIFC :=  mpif90 
  SCC :=  pgcc 
  SCXX :=  pgc++ 
  SFC :=  pgf95 
endif
NETCDF_PATH := $(NETCDF)
PIO_FILESYSTEM_HINTS := gpfs
PNETCDF_PATH := $(PNETCDF)
ifeq ($(COMPILER),intel)
  HAS_F2008_CONTIGUOUS := TRUE
  ifeq ($(MPILIB),mpi-serial)
    ifeq ($(compile_threaded),FALSE)
      PFUNIT_PATH := $(CESMDATAROOT)/tools/pFUnit/pFUnit3.2.8_cheyenne_Intel17.0.1_noMPI_noOpenMP
    endif
  endif
  ifeq ($(MPILIB),mpt)
    ifeq ($(compile_threaded),TRUE)
      PFUNIT_PATH := $(CESMDATAROOT)/tools/pFUnit/pFUnit3.2.8_cheyenne_Intel17.0.1_MPI_openMP
    endif
  endif
endif
CPPDEFS := $(CPPDEFS)  -DCESMCOUPLED 
ifeq ($(MODEL),pop)
  CPPDEFS := $(CPPDEFS)  -D_USE_FLOW_CONTROL 
endif
ifeq ($(MODEL),ufsatm)
  CPPDEFS := $(CPPDEFS)  -DSPMD 
  FFLAGS := $(FFLAGS)  $(FC_AUTO_R8) 
endif
ifeq ($(MODEL),gptl)
  CPPDEFS := $(CPPDEFS)  -DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY 
endif
ifeq ($(MODEL),mom)
  FFLAGS := $(FFLAGS)  $(FC_AUTO_R8) -Duse_LARGEFILE
endif
ifeq ($(COMPILER),gnu)
  FFLAGS := $(FFLAGS)  -fallow-argument-mismatch -fallow-invalid-boz 
  CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DNO_R16 -DCPRGNU
  SLIBS := $(SLIBS)  -ldl 
  ifeq ($(compile_threaded),TRUE)
    FFLAGS := $(FFLAGS)  -fopenmp 
    CFLAGS := $(CFLAGS)  -fopenmp 
  endif
  ifeq ($(DEBUG),TRUE)
    FFLAGS := $(FFLAGS)  -g -Wall -Og -fbacktrace -ffpe-trap=zero,overflow -fcheck=bounds 
    CFLAGS := $(CFLAGS)  -g -Wall -Og -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds 
  endif
  ifeq ($(DEBUG),FALSE)
    FFLAGS := $(FFLAGS)  -O 
    CFLAGS := $(CFLAGS)  -O 
  endif
  ifeq ($(MODEL),pio1)
    CPPDEFS := $(CPPDEFS)  -DNO_MPIMOD 
  endif
  ifeq ($(compile_threaded),TRUE)
    LDFLAGS := $(LDFLAGS)  -fopenmp 
  endif
endif
ifeq ($(COMPILER),intel)
  FFLAGS := $(FFLAGS)  -qopt-report -xCORE_AVX2 -no-fma
  CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DCPRINTEL
  CFLAGS := $(CFLAGS)  -qopt-report -xCORE_AVX2 -no-fma
  ifeq ($(compile_threaded),TRUE)
    FFLAGS := $(FFLAGS)  -qopenmp 
    CFLAGS := $(CFLAGS)  -qopenmp 
  endif
  ifeq ($(DEBUG),TRUE)
    FFLAGS := $(FFLAGS)  -O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created 
    CFLAGS := $(CFLAGS)  -O0 -g 
    CMAKE_OPTS := $(CMAKE_OPTS)  -DPIO_ENABLE_LOGGING=ON 
  endif
  ifeq ($(DEBUG),FALSE)
    FFLAGS := $(FFLAGS)  -O2 -debug minimal 
    CFLAGS := $(CFLAGS)  -O2 -debug minimal 
  endif
  ifeq ($(MPILIB),mpich)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpich2)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mvapich)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mvapich2)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpt)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),openmpi)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),impi)
    SLIBS := $(SLIBS)  -mkl=cluster 
  endif
  ifeq ($(MPILIB),mpi-serial)
    SLIBS := $(SLIBS)  -mkl 
  endif
  ifeq ($(compile_threaded),TRUE)
    LDFLAGS := $(LDFLAGS)  -qopenmp 
  endif
endif
ifeq ($(COMPILER),pgi)
  CPPDEFS := $(CPPDEFS)  -DFORTRANUNDERSCORE -DNO_SHR_VMATH -DNO_R16  -DCPRPGI 
  SLIBS := $(SLIBS)  -llapack -lblas 
  ifeq ($(compile_threaded),TRUE)
    FFLAGS := $(FFLAGS)  -mp 
  endif
  ifeq ($(MODEL),datm)
    FFLAGS := $(FFLAGS)  -Mnovect 
  endif
  ifeq ($(MODEL),dlnd)
    FFLAGS := $(FFLAGS)  -Mnovect 
  endif
  ifeq ($(MODEL),drof)
    FFLAGS := $(FFLAGS)  -Mnovect 
  endif
  ifeq ($(MODEL),dwav)
    FFLAGS := $(FFLAGS)  -Mnovect 
  endif
  ifeq ($(MODEL),dice)
    FFLAGS := $(FFLAGS)  -Mnovect 
  endif
  ifeq ($(MODEL),docn)
    FFLAGS := $(FFLAGS)  -Mnovect 
  endif
  ifeq ($(DEBUG),TRUE)
    FFLAGS := $(FFLAGS)  -O0 -g -Ktrap=fp -Mbounds -Kieee 
  endif
  ifeq ($(MPILIB),mpi-serial)
    SLIBS := $(SLIBS)  -ldl 
  endif
  ifeq ($(compile_threaded),TRUE)
    CFLAGS := $(CFLAGS)  -mp 
    LDFLAGS := $(LDFLAGS)  -mp 
  endif
endif
ifeq ($(MODEL),ufsatm)
  INCLDIR := $(INCLDIR)  -I$(EXEROOT)/atm/obj/FMS 
endif
