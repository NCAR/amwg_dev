#!/bin/bash

set -e

# Created 2024-06-07 15:13:13

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/b.e23_alpha17f.BLT1850.ne30_t232.098"

/glade/work/hannay/cesm_tags/cesm2_3_alpha17f/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./case.setup

./preview_namelists

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./preview_namelists

./case.build

./case.setup --reset

./preview_namelists

./preview_namelists

./case.build --clean-all

./preview_namelists

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./case.submit

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./preview_namelists

./case.build

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./case.build

./case.submit

./xmlchange MOM6_VERTICAL_GRID=hycom1

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./xmlchange RESUBMIT=10

./case.submit

./xmlchange RESUBMIT=15

./xmlchange RESUBMIT=10

./case.submit

./case.submit

./xmlchange RESUBMIT=10

./case.submit

./xmlchange RESUBMIT=10

./case.submit

./xmlchange RESUBMIT=10

./xmlchange RESUBMIT=50

./case.submit

