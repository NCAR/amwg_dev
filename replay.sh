#!/bin/bash

set -e

# Created 2023-07-27 14:31:26

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_119.FLTHIST_ne30.r328_gamma0.33_soae.001"

/glade/work/hannay/cesm_tags/cam6_3_119/cime/scripts/create_newcase --compset FLTHIST_v0d --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.setup

./xmlchange RUN_STARTDATE=1995-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1.001

./xmlchange RUN_REFDATE=1994-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./case.build

./case.build

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./xmlchange PROJECT=CESM0023,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths

./xmlchange REST_OPTION=nmonths,REST_N=1

./case.submit

./xmlchange PROJECT=CESM0023,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths

./xmlchange REST_OPTION=nmonths,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./xmlchange PROJECT=CESM0023,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./xmlchange RESUBMIT=2

./case.submit

