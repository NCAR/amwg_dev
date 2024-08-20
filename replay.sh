#!/bin/bash

set -e

# Created 2024-08-13 14:53:29

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/b.e30_beta02.BLTHIST.ne30_t232.104"

/glade/work/hannay/cesm_tags/cesm3_0_beta02/cime/scripts/create_newcase --compset BLTHIST --res ne30pg3_t232 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./case.setup

./preview_namelists

./xmlchange RUN_REFCASE=b.e23_alpha17f.BLT1850.ne30_t232.0104

./xmlchange RUN_REFDATE=0051-01-01

./xmlchange RUN_TYPE=hybrid

./xmlchange GET_REFCASE=true

./xmlchange RUN_REFDIR=cesm2_init

./xmlchange MOM6_VERTICAL_GRID=hycom1

./preview_namelists

./preview_namelists

./case.build

./case.build

./preview_namelists

./preview_namelists

./case.build

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./preview_namelists

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./xmlchange MOM6_VERTICAL_GRID=hycom1

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./xmlchange RUN_REFCASE=b.e30_beta02.BLT1850.ne30_t232.0104

./xmlchange RUN_REFDATE=0051-01-01

./xmlchange RUN_TYPE=hybrid

./xmlchange GET_REFCASE=true

./xmlchange RUN_REFDIR=cesm2_init

./preview_namelists

./preview_namelists

./case.build

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./preview_namelists

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./case.build

./case.submit

./case.build

./case.build

./case.submit

./xmlchange RUN_REFCASE=b.e30_beta02.BLT1850.ne30_t232.104

./xmlchange RUN_REFDATE=0051-01-01

./xmlchange RUN_TYPE=hybrid

./xmlchange GET_REFCASE=true

./xmlchange RUN_REFDIR=cesm2_init

./case.build

./case.submit

./case.build

./case.submit

./case.submit

./xmlchange RESUBMIT=20

