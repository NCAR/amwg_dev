#!/bin/bash

set -e

# Created 2022-01-25 14:32:39

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.e21.FWscHIST.ne30_L48_BL10_cam6_3_041_kzz3_zmtop50.hf.001"

/glade/work/hannay/cesm_tags/cam6_3_041/cime/scripts/create_newcase --compset FWscHIST --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -age_of_air_trcs -chem waccm_sc_mam4 -nlev 58

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

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./case.build

./case.build

./case.build --clean-all

./case.build

./case.build

./case.submit

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./case.submit

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./case.submit

./case.build

./xmlchange CONTINUE_RUN=TRUE

./xmlchange CONTINUE_RUN=FALSE

./case.submit

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./case.build

./case.submit

./xmlchange CONTINUE_RUN=FALSE

./case.submit

./case.submit

