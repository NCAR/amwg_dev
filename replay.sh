#!/bin/bash

set -e

# Created 2024-10-03 11:16:46

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_036.FLTHIST.ne120_L58.001"

/glade/work/hannay/cesm_tags/cam6_4_036/cime/scripts/create_newcase --compset FLTHIST --res ne120pg3_ne120pg3_mt13 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./xmlchange NTASKS=6144

./case.setup

./preview_namelists

./xmlchange RUN_STARTDATE=2000-01-01

./preview_namelists

./preview_namelists

./case.build

./preview_namelists

./case.build

./preview_namelists

./case.build

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=59,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=59,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_PRIORITY=premium

./case.submit

./case.submit

./xmlchange JOB_PRIORITY=regular

./case.submit

./xmlchange JOB_PRIORITY=premium

./case.submit

./case.submit

./case.submit

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./xmlchange NTASKS=6144

./case.setup

./preview_namelists

./xmlchange RUN_STARTDATE=2000-01-01

./case.submit

./xmlchange CAM_CONFIG_OPTS="-phys cam7 -nlev 58 -model_top lt -rad rrtmgp"

./case.build

./case.build

./case.build

./case.build

./case.submit

./case.submit

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./case.submit

./preview_namelists

./case.submit

./preview_namelists

./case.submit

./case.submit

./xmlchange CONTINUE_RUN=FALSE

./preview_namelists

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,RESUBMIT=59,STOP_N=2,STOP_OPTION=nmonths

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nmonths,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=premium

./case.submit

./case.build

