#!/bin/bash

set -e

# Created 2024-05-28 17:16:40

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/b.e23_alpha17f.BLT1850.ne30_t232.095"

/glade/work/hannay/cesm_tags/cesm2_3_alpha17f/cime/scripts/create_newcase --compset BLT1850_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./case.setup

./preview_namelists

./xmlchange MOM6_VERTICAL_GRID=hycom1

./preview_namelists

./preview_namelists

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./case.build

./case.build

./case.build

./case.build

./case.submit

./case.submit

./case.submit

./case.submit

./xmlchange RESUBMIT=20

./case.submit

./case.submit

./case.build

./case.submit

