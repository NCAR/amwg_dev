#!/bin/bash

set -e

# Created 2024-10-22 10:21:59

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_036.FLTHIST.ne30_L32.cam5.002"

/glade/work/hannay/cesm_tags/cam6_4_036/cime/scripts/create_newcase --compset FLTHIST --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS="-phys cam7 -nlev 32 -model_top lt -rad rrtmgp"

./xmlchange CAM_CONFIG_OPTS="-phys cam5 -nlev 32"

./xmlchange NTASKS=1280

./case.setup

./xmlchange RUN_STARTDATE=2000-01-01

./preview_namelists

./case.build

./case.submit

././xmlchange PROJECT=CESM0023,RESUBMIT=4,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=regular

./case.submit

./xmlchange CAM_CONFIG_OPTS="-phys cam5 -nlev 32"

./xmlchange NTASKS=1280

./case.setup

./xmlchange RUN_STARTDATE=2000-01-01

./preview_namelists

./preview_namelists

./preview_namelists

./case.build

./preview_namelists

./preview_namelists

./case.build

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=4,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=premium

./case.submit

