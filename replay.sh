#!/bin/bash

set -e

# Created 2022-06-10 12:07:54

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cesm3_cam058_mom.FWscHIST.ne30_L58.noclubbtop-plus.001"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_058_MOM/cime/scripts/create_newcase --compset FWscHIST --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722 --driver nuopc

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2 -chem waccm_sc_mam4 -nlev 58

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

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=3,STOP_N=3,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

