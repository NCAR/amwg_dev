#!/bin/bash

set -e

# Created 2024-05-31 11:37:00

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/b.e23_alpha17f.BLTHIST.ne30_t232.092"

/glade/work/hannay/cesm_tags/cesm2_3_alpha17f/cime/scripts/create_newcase --compset BLTHIST_v0c --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./case.setup

./preview_namelists

./xmlchange MOM6_VERTICAL_GRID=hycom1

./preview_namelists

./case.setup

./preview_namelists

./xmlchange RUN_REFCASE=b.e23_alpha17f.BLTHIST.ne30_t232.092

./xmlchange RUN_REFDATE=0101-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./preview_namelists

./case.build

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=b.e23_alpha17f.BLTHIST.ne30_t232.092

./xmlchange RUN_REFDATE=0101-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./case.build

./case.build

./case.build

./xmlchange MOM6_VERTICAL_GRID=hycom1

./preview_namelists

./case.setup

./preview_namelists

./xmlchange RUN_TYPE=hybrid

./xmlchange RUN_REFCASE=b.e23_alpha17f.BLT1850.ne30_t232.092

./xmlchange RUN_REFDATE=0101-01-01

./xmlchange GET_REFCASE=TRUE

./xmlchange RUN_REFDIR=cesm2_init

./xmlchange MOM6_VERTICAL_GRID=hycom1

./preview_namelists

./preview_namelists

./preview_namelists

./preview_namelists

./preview_namelists

./case.build

./preview_namelists

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=6:00:00

./xmlchange USER_REQUESTED_WALLTIME=12:00:00

./case.submit

./case.submit

./case.submit

./case.submit

./case.submit

./preview_namelists

./preview_namelists

./case.submit

./preview_namelists

./case.submit

./xmlchange RESUBMIT=20

./xmlchange RESUBMIT=40

./case.submit

./case.submit

./xmlchange CONTINUE_RUN=25

./case.submit

./xmlchange RESUBMIT=15

./xmlchange RESUBMIT=15

