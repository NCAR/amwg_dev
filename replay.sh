#!/bin/bash

set -e

# Created 2024-10-16 10:37:29

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_036.FLTHIST.ne30_L58.SST4K.001"

/glade/work/hannay/cesm_tags/cam6_4_036/cime/scripts/create_newcase --compset FLTHIST --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./xmlchange NTASKS=2160

./case.setup

./xmlchange RUN_STARTDATE=2000-01-01

./xmlchange SSTICE_DATA_FILENAME=/glade/campaign/cesm/cesmdata/cseg/inputdata/atm/cam/sst/sst_HadOIBl_bc_1x1_1850_2021_SST4K_c241009.nc

./preview_namelists

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=4,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=premium

./case.submit

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=4,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=premium

./case.submit

./case.submit

./xmlchange JOB_PRIORITY=regular

