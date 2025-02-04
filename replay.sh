#!/bin/bash

set -e

# Created 2025-01-12 10:23:21

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_beta04.FLTHIST.ne30.hetfrz_dust.001"

/glade/work/hannay/cesm_tags/cesm3_0_beta04/cime/scripts/create_newcase --compset FLTHIST --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS="-rad rrtmgp"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange NTASKS=2160

./case.setup

./xmlchange RUN_STARTDATE=2000-01-01

./preview_namelists

./preview_namelists

./case.build

./preview_namelists

./case.build

./preview_namelists

./case.build

./preview_namelists

./case.build

./preview_namelists

./case.build

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=4,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=2

./xmlchange JOB_WALLCLOCK_TIME=12:00:00 --subgroup case.run

./xmlchange JOB_WALLCLOCK_TIME=06:00:00 --subgroup case.st_archive

./xmlchange JOB_PRIORITY=regular

./case.submit

