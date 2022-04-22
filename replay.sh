#!/bin/bash

set -e

# Created 2022-04-22 15:12:51

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.cesm3_cam041_mom.B1850MOM.ne30_L58_t061.003"

/glade/work/hannay/cesm_tags/cesm3_cam6_3_041_MOM3/cime/scripts/create_newcase --compset B1850MOM --res ne30pg3_t061 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange CAM_CONFIG_OPTS=-phys cam_dev -microphys mg2 -chem waccm_sc_mam4 -nlev 58

./xmlchange ATM_GRID=ne30np4.pg3

./case.setup

./case.build

./case.build

./case.build

./case.build

./case.build

./case.build

./xmlchange PROJECT=P93300642,JOB_QUEUE=regular,RESUBMIT=1,STOP_N=1,STOP_OPTION=nmonths

./xmlchange JOB_QUEUE=premium

./case.submit

./case.build

./case.build

./case.build

./case.build

./case.submit

./xmlchange PROJECT=CESM0019,JOB_QUEUE=regular,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./case.submit

./case.setup --reset

./case.build

./case.build --clean-all

./case.build

