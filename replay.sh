#!/bin/bash

set -e

# Created 2023-08-29 09:42:33

CASEDIR="/glade/p/cesmdata/cseg/runs/cesm2_0/b.e23_alpha16b.BLT1850.ne30_t232.040"

/glade/work/hannay/cesm_tags/cesm2_3_alpha16b/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./case.build

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./preview_namelists

./case.build

./case.build --clean-all

./case.build

./case.build

./xmlchange RUN_REFCASE=b.e23_alpha16b.BLT1850.ne30_t232.040

./xmlchange RUN_REFDATE=0021-01-01

./xmlchange RUN_STARTDATE=0021-01-01

./xmlchange RUN_TYPE=hybrid

./case.build

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=0,STOP_N=1,STOP_OPTION=nmonths

./xmlchange REST_OPTION=nmonths,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

./xmlchange RUN_REFCASE=b.e23_alpha16b.BLT1850.ne30_t232.040

./xmlchange RUN_REFDATE=0021-01-01

./xmlchange RUN_STARTDATE=0021-01-01

./xmlchange RUN_TYPE=branch

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00

./case.submit

