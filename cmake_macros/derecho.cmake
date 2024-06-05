if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_NANOTIME -DBIT64 -DHAVE_VPRINTF -DHAVE_BACKTRACE -DHAVE_SLASHPROC -DHAVE_COMM_F2C -DHAVE_TIMES -DHAVE_GETTIMEOFDAY")
endif()
set(MPI_SERIAL_PATH "$ENV{NCAR_ROOT_MPI_SERIAL}")
set(NETCDF_PATH "$ENV{NETCDF}")
set(PIO_FILESYSTEM_HINTS "lustre")
set(PNETCDF_PATH "$ENV{PNETCDF}")
if(DEFINED ENV{PIO})
  set(PIO_LIBDIR "$ENV{PIO}/lib")
  set(PIO_INCDIR "$ENV{PIO}/include")
endif()
# If we want to use cray-libsci instead of mkl uncomment this line as well as the module in config_machines.xml
string(REPLACE "-mkl=cluster" "" SLIBS "${SLIBS}")
string(APPEND CPPDEFS " -DHAVE_GETTID")
