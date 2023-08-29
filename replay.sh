#!/bin/bash

set -e

# Created 2023-08-25 16:26:12

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.e23_alpha16b.BLT1850.ne30_t232.041"

/glade/work/hannay/cesm_tags/cesm2_3_alpha16b/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.build

./case.build

./preview_namelists

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./xmlchange RESUBMIT=0

./xmlchange RESUBMIT=0

