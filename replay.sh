#!/bin/bash

set -e

# Created 2023-04-24 16:47:39

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1_ebudget_tuning1.001"

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

./xmlchange RUN_STARTDATE=1995-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1_ebudget.001

./xmlchange RUN_REFDATE=1995-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./case.build

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

./xmlchange RUN_STARTDATE=1995-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=f.cam6_3_107.FLTHIST_v0a.ne30.clm5_1.001

./xmlchange RUN_REFDATE=1986-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./preview_namelists

./preview_namelists

./case.build

./case.build

./case.build --clean-all

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./xmlchange JOB_QUEUE=premium

./case.submit

./xmlchange JOB_QUEUE=regular

