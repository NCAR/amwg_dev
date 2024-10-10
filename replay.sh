#!/bin/bash

set -e

# Created 2024-09-06 14:54:20

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_beta02.FLTHIST.ne30.103"

/glade/work/hannay/cesm_tags/cesm3_0_beta02/cime/scripts/create_newcase --compset FLTHIST --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722 --pecount 2160

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./case.setup

./preview_namelists

./xmlchange RUN_REFCASE=b.e30_beta02.BLTHIST.ne30_t232.104

./xmlchange RUN_REFDATE=1900-01-01

./xmlchange RUN_TYPE=hybrid

./xmlchange GET_REFCASE=true

./xmlchange RUN_REFDIR=cesm2_init

./xmlchange RUN_STARTDATE=1900-01-01

./preview_namelists

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=10,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./case.submit

./case.build

./case.submit

./case.submit

./case.build

./case.submit

./xmlchange RUN_TYPE=startup

./preview_namelists

./case.submit

./case.submit

./case.submit

./preview_namelists

./case.submit

./case.submit

./preview_namelists

./case.submit

./case.submit

./xmlchange RESUBMIT=10

./xmlchange RESUBMIT=22

./case.submit

