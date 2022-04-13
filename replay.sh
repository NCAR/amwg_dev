#!/bin/bash

set -e

# Created 2022-04-12 16:36:22

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.c6_3_41.FWscHIST.ne30_L58.zm2_fix.002"

/glade/work/hannay/cesm_tags/cam6_3_041/cime/scripts/create_newcase --compset FWscHIST --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev  -chem waccm_sc_mam4 -nlev 58

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

./xmlchange JOB_QUEUE=premium,STOP_N=1,STOP_OPTION=nmonths,RESUBMIT=0

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=0,STOP_N=1,STOP_OPTION=nyears

./case.submit

