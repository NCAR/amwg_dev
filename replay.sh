#!/bin/bash

set -e

# Created 2023-04-20 14:04:12

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1_ebudget.001"

/glade/work/hannay/cesm_tags/cam6_3_107_ebudget_dev_update/cime/scripts/create_newcase --compset HIST_CAM%DEV%LT%GHGMAM4_CLM51%SP_CICE%PRES_DOCN%DOM_MOSART_SGLC_SWAV --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

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

./preview_namelists

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=premium,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./preview_namelists

./preview_namelists

./xmlchange PROJECT=P93300642,JOB_QUEUE=premium,RESUBMIT=1,STOP_N=1,STOP_OPTION=ndays

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./case.submit

./xmlchange PROJECT=P93300642,JOB_QUEUE=premium,RESUBMIT=1,STOP_N=1,STOP_OPTION=ndays

./xmlchange REST_OPTION=ndays,REST_N=1

./case.submit

./xmlchange CONTINUE_RUN=FALSE

./case.submit

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=15,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

