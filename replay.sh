#!/bin/bash

set -e

# Created 2024-06-04 21:07:47

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/b.e23_alpha17f.BLT1850.ne30_t232.097"

/glade/work/hannay/cesm_tags/cesm2_3_alpha17f/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange MOM6_VERTICAL_GRID=hycom1

./preview_namelists

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./preview_namelists

./case.submit

./case.build

./case.submit

./case.submit

