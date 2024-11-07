#!/bin/bash

set -e

# Created 2024-11-06 14:17:53

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_036.FLTHIST.ne30_L32.cam5.rpen.003"

/glade/work/hannay/cesm_tags/cam6_4_036/cime/scripts/create_newcase --compset FLTHIST --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS="-phys cam5 -nlev 32"

./xmlchange NTASKS=1280

./case.setup

./xmlchange RUN_STARTDATE=2000-01-01

./preview_namelists

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=4,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=regular

./case.submit

