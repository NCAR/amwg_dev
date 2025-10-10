#!/bin/bash

set -e

# Created 2025-10-10 15:53:34

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_116.FHISTC_LTso.ne120_L58.test002"

/glade/work/hannay/cesm_tags/cam6_4_116/cime/scripts/create_newcase --compset FHISTC_LTso --res ne120pg3_ne120pg3_mt13 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange NTASKS=6144

./case.setup

./preview_namelists

./xmlchange CASE_GIT_REPOSITORY=git@github.com:NCAR/amwg_dev.git

./xmlchange RUN_STARTDATE=1990-01-01

./preview_namelists

./preview_namelists

./case.build

./preview_namelists

./case.build

./xmlchange SSTICE_DATA_FILENAME=$DIN_LOC_ROOT/atm/cam/sst/sst_HadOIBl_bc_1x1_1850_2021_c120422.nc

./xmlchange SSTICE_DATA_FILENAME=/glade/work/xueliu/data/cesm_input/SST/SSTinput_CDR_noleap/sst_input_0.25_CDR-19900101-19901231.nc

./xmlchange SSTICE_MESH_FILENAME=/glade/campaign/cgd/amp/juliob/NOAA_OI_SST/sst_ice_NOAA_QxQ_ESMFmesh.nc

./xmlchange SSTICE_YEAR_START=1990

./xmlchange SSTICE_YEAR_END=1990

./xmlchange SSTICE_YEAR_ALIGN=1990

./preview_namelists

./case.build

./case.build --clean-all

