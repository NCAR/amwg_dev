#!/bin/bash

set -e

# Created 2024-03-07 14:21:55

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/b.e23_alpha16g.BLT1850.ne30_t232.078d"

/glade/work/hannay/cesm_tags/cesm2_3_alpha16g/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange RUN_REFCASE=b.e23_alpha16g.BLT1850.ne30_t232.075

./xmlchange RUN_REFDATE=0075-01-01

./xmlchange RUN_TYPE=hybrid

./xmlchange GET_REFCASE=true

./xmlchange RUN_REFDIR=cesm2_init

./preview_namelists

./preview_namelists

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=24,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./case.setup

./case.build --clean

./preview_namelists

./xmlchange RUN_REFCASE=b.e23_alpha16g.BLT1850.ne30_t232.075

./xmlchange RUN_REFDATE=0075-01-01

./xmlchange RUN_TYPE=hybrid

./xmlchange GET_REFCASE=true

./xmlchange RUN_REFDIR=cesm2_init

./preview_namelists

./case.build

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./xmlchange PROJECT=CESM0023,RESUBMIT=24,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=24,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=24,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./case.submit

./xmlchange RUN_REFCASE=b.e23_alpha16g.BLT1850.ne30_t232.075

./xmlchange RUN_REFDATE=0075-01-01

./xmlchange RUN_TYPE=hybrid

./xmlchange GET_REFCASE=true

./xmlchange RUN_REFDIR=cesm2_init

./case.submit

./case.submit

./xmlchange RESUBMIT=0

