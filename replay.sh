#!/bin/bash

set -e

# Created 2024-08-22 16:47:57

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/b.e30_beta02.BLT1850.ne30_t232.109"

/glade/work/hannay/cesm_tags/cesm3_0_beta02/cime/scripts/create_newcase --compset BLT1850 --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722 --pecount L

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./case.setup

./preview_namelists

./preview_namelists

./case.build

./xmlchange RUN_REFCASE=b.e23_alpha17f.BLT1850.ne30_t232.098

./xmlchange RUN_REFDATE=0201-01-01

./xmlchange RUN_TYPE=hybrid

./xmlchange GET_REFCASE=true

./xmlchange RUN_REFDIR=cesm2_init

./xmlchange MOM6_VERTICAL_GRID=hycom1

./preview_namelists

./case.build

./preview_namelists

./xmlchange MOM6_VERTICAL_GRID=hycom1

./preview_namelists

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

./case.build --clean-all

./preview_namelists

./case.build

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./xmlchange RESUBMIT=20

./xmlchange RESUBMIT=25

./case.submit

