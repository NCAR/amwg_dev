#!/bin/bash

set -e

# Created 2024-01-09 11:49:32

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.e23_alpha16g.BLT1850.ne30_t232.064intel"

/glade/work/hannay/cesm_tags/cesm2_3_alpha16g/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722 --compiler intel-oneapi

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./preview_namelists

./case.build

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=12,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./xmlchange RESUBMIT=20

./xmlchange RESUBMIT=0

./case.submit

./xmlchange RESUBMIT=12

