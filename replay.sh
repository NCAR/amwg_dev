#!/bin/bash

set -e

# Created 2025-08-15 17:30:08

CASEDIR="/glade/campaign/cesm/cesmdata/cseg/runs/cesm2_0/f.e30_cam6_4_097.FHISTC_MTso.ne30_greendlndantarcsgh30fac2.5_nu_div.001"

/glade/work/hannay/cesm_tags/cam6_4_097/cime/scripts/create_newcase --compset FHISTC_MTso --res ne30pg3_ne30pg3_mg17 --case "${CASEDIR}" --run-unsupported --project 93300722

cd "${CASEDIR}"

./xmlchange --append CAM_CONFIG_OPTS=-cosp

./xmlchange NTASKS=2160

./case.setup

./xmlchange RUN_STARTDATE=1990-01-01

./preview_namelists

./xmlchange CASE_GIT_REPOSITORY=git@github.com:NCAR/amwg_dev.git

./preview_namelists

./case.build

./xmlchange PROJECT=CESM0023,RESUBMIT=6,STOP_N=2,STOP_OPTION=nyears

./xmlchange CHARGE_ACCOUNT=CESM0023,PROJECT=CESM0023

./xmlchange REST_OPTION=nyears,REST_N=1

