#!/bin/bash

set -e

# Created 2025-10-11 10:39:17

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_116.FHISTC_LTso.ne120_L58.test003"

/glade/work/hannay/cesm_tags/cam6_4_116/cime/scripts/create_newcase --compset FHISTC_LTso --res ne120pg3_ne120pg3_mt13 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange NTASKS=6144

./case.setup

./preview_namelists

./xmlchange CASE_GIT_REPOSITORY=git@github.com:NCAR/amwg_dev.git

./xmlchange RUN_STARTDATE=1990-01-01

./xmlchange RUN_STARTDATE=1990-01-01

./preview_namelists

./preview_namelists

./case.build

./xmlchange SSTICE_DATA_FILENAME=glade/derecho/scratch/juliob/SSTinput_CDR_noleap/decades_cdf5/sst_input_0.25_CDR-19890101-20011231.nc

./xmlchange SSTICE_MESH_FILENAME=/glade/campaign/cgd/amp/juliob/NOAA_OI_SST/sst_ice_NOAA_QxQ_ESMFmesh.nc

./xmlchange SSTICE_YEAR_START=1989

./xmlchange SSTICE_YEAR_END=2001

./xmlchange SSTICE_YEAR_ALIGN=1850

./xmlchange DOUT_S_ROOT=/glade/campaign/cesm/community/amwg/cam_hr/$CASE

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=1,STOP_N=1,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=premium

./xmlchange JOB_PRIORITY=special --force

./case.submit

./xmlchange SSTICE_DATA_FILENAME=glade/derecho/scratch/juliob/SSTinput_CDR_noleap/decades_cdf5/sst_input_0.25_CDR-19890101-20011231.nc

./xmlchange SSTICE_MESH_FILENAME=/glade/campaign/cgd/amp/juliob/NOAA_OI_SST/sst_ice_NOAA_QxQ_ESMFmesh.nc

./xmlchange SSTICE_YEAR_START=1989

./xmlchange SSTICE_YEAR_END=2001

./xmlchange SSTICE_YEAR_ALIGN=1989

./xmlchange DOUT_S_ROOT=/glade/campaign/cesm/community/amwg/cam_hr/$CASE

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=1,STOP_N=1,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=premium

./xmlchange JOB_PRIORITY=special --force

./case.submit

./xmlchange DOUT_S_ROOT=/glade/campaign/cesm/community/amwg/cam_hr/$CASE

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=1,STOP_N=1,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=premium

./xmlchange JOB_PRIORITY=special --force

./case.submit

./xmlchange SSTICE_DATA_FILENAME=/glade/derecho/scratch/juliob/SSTinput_CDR_noleap/decades_cdf5/sst_input_0.25_CDR-19890101-20011231.nc

