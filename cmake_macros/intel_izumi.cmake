if (MPILIB STREQUAL mpi-serial AND NOT compile_threaded)
  set(PFUNIT_PATH "/fs/cgd/csm/tools/pFUnit/pFUnit4.7.0_izumi_Intel20.0.1_noMPI_noOpenMP")
endif()
