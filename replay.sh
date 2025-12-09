#!/bin/bash

set -e

# Created 2025-12-09 14:55:24

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_134.FHISTC_LTso.ne30_L58.001"

/glade/work/hannay/cesm_tags/cam6_4_134/cime/scripts/create_newcase --compset FHISTC_LTso --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange NTASKS=2176

./case.setup

./preview_namelists

./xmlchange CASE_GIT_REPOSITORY=git@github.com:NCAR/amwg_dev.git

./xmlchange RUN_STARTDATE=1990-01-01

./preview_namelists

./case.build

./xmlchange NTASKS=2700

./case.setup --reset

./xmlchange DOUT_S_ROOT=/glade/campaign/cesm/community/amwg/cam_hr/$CASE

./xmlchange PROJECT=CESM0023,RESUBMIT=9,STOP_N=4,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1,DOUT_S_SAVE_INTERIM_RESTART_FILES=TRUE

./xmlchange JOB_WALLCLOCK_TIME=24:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=premium

./preview_namelists

./case.build

./xmlchange DOUT_S_ROOT=/glade/campaign/cesm/community/amwg/cam_hr/$CASE

