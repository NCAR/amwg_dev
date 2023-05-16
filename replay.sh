#!/bin/bash

set -e

# Created 2023-05-10 15:15:14

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_110.FMTHIST_v0c.ne30.tuningF.001"

/glade/work/hannay/cesm_tags/cam6_3_110/cime/scripts/create_newcase --compset FMTHIST_v0c --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

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

./xmlchange PROJECT=P93300642,JOB_QUEUE=premium,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./case.submit

./xmlchange STOP_N=nyears,RESUBMIT=10

./xmlchange STOP_N=1,STOP_OPTION=nyears,RESUBMIT=10

./case.submit

./xmlchange STOP_N=1,STOP_OPTION=nyears,RESUBMIT=10,JOB_QUEUE=regular

./case.submit

