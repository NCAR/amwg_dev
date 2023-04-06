#!/bin/bash

set -e

# Created 2023-04-03 13:44:32

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_101.FLTHIST_v0a.ne30.003"

/glade/work/hannay/cesm_tags/cam6_3_101/cime/scripts/create_newcase --compset FLTHIST_v0a --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.setup

./xmlchange RUN_STARTDATE=1979-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=f.e21.FWscHIST_BGC.ne30_ne30_mg17_L48_revert-J.001

./xmlchange RUN_REFDATE=1989-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./case.build

./case.build

./case.build

./case.build

./case.build --clean-all

./case.build

./case.build

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=premium,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

