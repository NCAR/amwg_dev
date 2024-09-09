#!/bin/bash

set -e

# Created 2024-04-24 21:15:59

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/b.e23_alpha16g.BLT1850.ne30_t232_daily.082b"

/glade/work/hannay/cesm_tags/cesm2_3_alpha16g/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange RUN_REFCASE=b.e23_alpha16g.BLT1850.ne30_t232.075

./xmlchange RUN_REFDATE=0101-01-01

./xmlchange RUN_TYPE=hybrid

./xmlchange GET_REFCASE=true

./xmlchange RUN_REFDIR=cesm2_init

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=1,STOP_N=1,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./preview_namelists

./case.build

./case.build

./case.submit

