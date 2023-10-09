#!/bin/bash

set -e

# Created 2023-09-29 11:16:53

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.e23_alpha16b.BLT1850.ne30_t232.049"

/glade/work/hannay/cesm_tags/cesm2_3_alpha16b_taus/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.build

./case.build

./case.build

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./xmlchange JOB_QUEUE=regular

./xmlchange JOB_QUEUE=premium

./case.submit

./xmlchange NTASKS_LND=720,NTASKS_ROF-720,ROOTPE_ICE=720

./xmlchange NTASKS_LND=720,NTASKS_ROF=720,ROOTPE_ICE=720

./case.setup --reset

./case.build --clean-all

./case.build

./case.submit

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./xmlchange JOB_QUEUE=regular

./case.submit

./case.submit

./case.submit

./case.submit

./xmlchange RESUBMIT=10

