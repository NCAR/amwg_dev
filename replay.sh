#!/bin/bash

set -e

# Created 2024-07-31 13:44:00

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e23_beta02.FLTHIST_ne30.001"

/glade/work/tilmes/my_cesm_sandbox/cesm3_0_beta02/cime/scripts/create_newcase --compset FLTHIST --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --pecount 2160 --project P93300643

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./xmlchange RUN_STARTDATE=1995-01-01

./xmlchange STOP_N=2

./xmlchange STOP_OPTION=nyears

./xmlchange RESUBMIT=5

./xmlchange TIMER_LEVEL=10

./xmlchange RUN_STARTDATE=1995-01-01

./case.setup

./case.setup

./case.build

./case.submit

./case.submit

./case.submit

./case.submit

