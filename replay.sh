#!/bin/bash

set -e

# Created 2023-05-26 13:09:11

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.e23_alpha14a.BLT1850.ne30_t232.032"

/glade/work/hannay/cesm_tags/cesm2_3_alpha14a/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange NTASKS=1800

./xmlchange NTASKS_OCN=324

./xmlchange NTASKS_WAV=36

./xmlchange NTASKS_GLC=36

./xmlchange NTASKS_ICE=1080

./xmlchange NTASKS_ROF=724

./xmlchange NTASKS_LND=724

./xmlchange NTASKS_ESP=1

./xmlchange ROOTPE_OCN=1800

./xmlchange ROOTPE_WAV=1

./xmlchange ROOTPE_ICE=724

./xmlchange NTASKS_LND=720,NTASKS_ROF=720,ROOTPE_ICE=720

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.build

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,JOB_QUEUE=premium,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./xmlchange PROJECT=CESM0023,JOB_QUEUE=premium,RESUBMIT=10,STOP_N=1,STOP_OPTION=nyears

./xmlchange PROJECT=CESM0023,JOB_QUEUE=regular,RESUBMIT=15,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./case.submit

./xmlchange STOP_N=1

./case.submit

./case.submit

./case.submit

./case.submit

./case.submit

