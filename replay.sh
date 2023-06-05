#!/bin/bash

set -e

# Created 2023-05-24 15:35:46

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/f.cam6_3_112.FLTHIST_v0c.ne30_L32.non-ogw-ubcT.001"

/glade/work/hannay/cesm_tags/cam6_3_112/cime/scripts/create_newcase --compset FLTHIST_v0c --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange --append CAM_CONFIG_OPTS=-phys cam_dev -chem ghg_mam4 -nlev 58 -cosp

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

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./xmlchange PROJECT=P93300642,JOB_QUEUE=premium,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -chem ghg_mam4 -nlev 32 -cosp

./case.setup --reset

./case.build --clean-all

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=premium,RESUBMIT=5,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./preview_namelists

./preview_namelists

./preview_namelists

./xmlchange JOB_QUEUE=regular

