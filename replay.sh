#!/bin/bash

set -e

# Created 2023-05-05 15:19:52

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_109.FLTHIST_v0b.ne30.tuningF.001"

/glade/work/hannay/cesm_tags/cam6_3_109/cime/scripts/create_newcase --compset FLTHIST_v0b --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.setup

./xmlchange RUN_STARTDATE=1995-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1.001

./xmlchange RUN_REFDATE=1995-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./preview_namelists

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

./preview_namelists

./preview_namelists

./preview_namelists

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=premium,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=premium,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./preview_namelists

./case.submit

./xmlchange STOP_N=2,STOP_OPTION=nyears,JOB_QUEUE=regular

./xmlchange STOP_N=2,STOP_OPTION=nyears,JOB_QUEUE=regular,RESUBMIT=5

./preview_namelists

./case.submit

./xmlchange STOP_N=1,STOP_OPTION-nmonths

./xmlchange STOP_N=1,STOP_OPTION=nmonths

./case.submit

./xmlchange RESUBMIT=10

./xmlchange CONTINUE_RUN=FALSE

./case.submit

./xmlchange STOP_OPTION=nyears

./xmlchange CONTINUE_RUN=FALSE

./case.submit

